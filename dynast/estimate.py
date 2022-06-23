import datetime as dt
import os
from typing import Dict, List, Optional, Union

import pandas as pd
import pystan
from typing_extensions import Literal

from . import config, constants, estimation, preprocessing, utils
from .logging import logger
from .stats import Stats


@logger.namespaced('estimate')
def estimate(
    count_dirs: List[str],
    out_dir: str,
    reads: Union[Literal['complete'], List[Literal['total', 'transcriptome', 'spliced', 'unspliced']]] = 'complete',
    groups: Optional[List[Dict[str, str]]] = None,
    ignore_groups_for_est: bool = True,
    genes: Optional[List[str]] = None,
    downsample: Optional[Union[int, float]] = None,
    downsample_mode: Literal['uniform', 'cell', 'group'] = 'uniform',
    cell_threshold: int = 1000,
    cell_gene_threshold: int = 16,
    control_p_e: Optional[float] = None,
    control: bool = False,
    method: Literal["pi_g", "alpha"] = "alpha",
    n_threads: int = 8,
    temp_dir: Optional[str] = None,
    nasc: bool = False,
    by_name: bool = False,
    seed: Optional[int] = None,
):
    """Main interface for the `estimate` command.

    Args:
        count_dirs: Paths to directories containing `count` command output
        out_dir: Output directory
        reads: What read group(s) to quantify
        groups: Cell groups
        ignore_groups_for_est: Ignore cell groups for final estimation
        genes: Genes to consider
        downsample: Downsample factor (float) or number (int)
        donsample_mode: Downsampling mode
        cell_threshold: Run estimation only for cells with at least this many counts
        cell_gene_threshold: Run estimation for cell-genes with at least this many counts
        control_p_e: Old RNA conversion rate (p_e), estimated from control samples
        control: Whether this is a control sample
        method: Estimation method to use
        n_threads: Number of threads
        temp_dir: Temporary directory
        nasc: Whether to match NASC-seq pipeline behavior
        by_name: Whether to group counts by gene name instead of ID
        seed: Random seed
    """
    stats = Stats()
    stats.start()
    stats_path = os.path.join(
        out_dir, f'{constants.STATS_PREFIX}_{dt.datetime.strftime(stats.start_time, "%Y%m%d_%H%M%S_%f")}.json'
    )
    os.makedirs(out_dir, exist_ok=True)

    # Check that all the conversions are the same if there are multiple count dirs
    conversions = utils.read_pickle(os.path.join(count_dirs[0], constants.CONVS_FILENAME))
    for count_dir in count_dirs[1:]:
        _conversions = utils.read_pickle(os.path.join(count_dir, constants.CONVS_FILENAME))
        if conversions != _conversions:
            raise Exception(
                f'Conversions for {count_dir} doesn\'t match conversions for {count_dirs[0]}. '
                f'({_conversions} != {conversions}).'
            )
    logger.info(f'Conversions: {" ".join(",".join(convs) for convs in conversions)}')
    all_conversions = sorted(utils.flatten_iter(conversions))

    gene_infos = utils.read_pickle(os.path.join(count_dirs[0], constants.GENES_FILENAME))

    # Read each counts dataframe and suffix barcodes if needed
    dfs = []
    for i, count_dir in enumerate(count_dirs):
        counts_path = os.path.join(count_dir, f'{constants.COUNTS_PREFIX}_{"_".join(all_conversions)}.csv')
        logger.info(
            f'Reading {counts_path}' + (f' and suffixing all barcodes with `-{i}`' if len(count_dirs) > 1 else '')
        )
        _df_counts = preprocessing.read_counts(counts_path)
        # Convert gene IDs to names
        if by_name:
            _df_counts['GX'] = _df_counts['GX'].apply(lambda gx: gene_infos[gx]['gene_name'] or gx)
        # Subset to provided genes
        if genes:
            _df_counts = _df_counts[_df_counts['GX'].isin(genes)]
        if len(count_dirs) > 1:
            _df_counts['barcode'] = _df_counts['barcode'].astype(str) + f'-{i}'
        dfs.append(_df_counts)
    df_counts_uncomplemented = pd.concat(dfs, ignore_index=True) if len(count_dirs) > 1 else dfs[0]
    df_counts_uncomplemented['barcode'] = df_counts_uncomplemented['barcode'].astype('category')
    df_counts_uncomplemented.drop(columns=['umi'], inplace=True)

    # Clean groups by combining multiple into one
    if isinstance(groups, list):
        if len(groups) == 1:
            groups = groups[0]
        else:
            groups = {f'{barcode}-{i}': group for i, _groups in enumerate(groups) for barcode, group in _groups.items()}
        # Contains group name to list of cells mapping
        group_cells = {}
        for barcode, group in groups.items():
            group_cells.setdefault(group, []).append(barcode)

    # Change barcodes to cell groups. We add groups here instead of after
    # complementing because the user may want to downsample per group instead of per cell.
    if groups:
        logger.warning(f'Barcodes that are not among the {len(groups)} barcodes with assigned groups will be ignored.')
        df_counts_uncomplemented = df_counts_uncomplemented[df_counts_uncomplemented['barcode'].isin(
            groups.keys()
        )].reset_index(drop=True)
        df_counts_uncomplemented['group'] = df_counts_uncomplemented['barcode'].map(groups).astype('category')

    # Downsample here
    if downsample:
        _proportion = None
        _count = None
        _group_by = None
        if int(downsample) == downsample:
            _count = int(downsample)
        else:
            _proportion = downsample
        if downsample_mode == 'cell':
            _group_by = ['barcode']
        elif downsample_mode == 'group':
            _group_by = ['group']
        if not _group_by:
            logger.info(
                'Downsampling uniformly at random ' +
                (f'to {_count} entries' if _count else f'to a factor of {_proportion}')
            )
        else:
            logger.info(f'Downsampling per {_group_by} to {_count} entries')
        old_count = df_counts_uncomplemented.shape[0]
        df_counts_uncomplemented = utils.downsample_counts(
            df_counts_uncomplemented, proportion=_proportion, count=_count, seed=seed, group_by=_group_by
        )
        logger.debug(f'Downsampled from {old_count} to {df_counts_uncomplemented.shape[0]} entries')

    # Check that all requested read groups can be corrected.
    transcriptome_any = df_counts_uncomplemented['transcriptome'].any()
    transcriptome_all = df_counts_uncomplemented['transcriptome'].all()
    if reads != 'complete':
        if 'transcriptome' in reads and not transcriptome_any:
            raise Exception(
                'No reads are assigned to `transcriptome`, so estimation is not supported for this read group.'
            )
        if 'total' in reads and transcriptome_all:
            raise Exception(
                'All reads are assigned to `transcriptome`, so estimation is not supported for `total` read group.'
            )
        for key in set(reads).intersection(('spliced', 'unspliced')):
            if not (df_counts_uncomplemented['velocity'] == key).any():
                raise Exception(
                    f'No reads are assigned to `{key}`, so estimation is not supported for this read group.'
                )
    else:
        reads = []
        if transcriptome_any:
            reads.append('transcriptome')
        if not transcriptome_all:
            reads.append('total')
        reads += list(set(df_counts_uncomplemented['velocity'].unique()) - set(config.VELOCITY_BLACKLIST))
    logger.info(f'Estimation will be done on the following read groups: {reads}')

    df_counts = preprocessing.complement_counts(df_counts_uncomplemented, gene_infos)
    logger.info(
        f'Final counts: {df_counts.shape[0]} reads '
        f'across {df_counts["barcode"].nunique()} barcodes' +
        (f' and {df_counts["group"].nunique()} groups.' if groups else '.')
    )

    # Estimate p_e
    p_key = 'group' if groups else 'barcode'
    p_e_path = os.path.join(out_dir, constants.P_E_FILENAME)
    if control_p_e:
        logger.info('`--p-e` provided. No background mutation rate estimation will be done.')
        df_barcodes = df_counts[[p_key]].drop_duplicates().reset_index(drop=True)
        df_barcodes['p_e'] = control_p_e
        df_barcodes.to_csv(p_e_path, header=[p_key, 'p_e'], index=False)
    else:
        logger.info(f'Estimating average conversion rate in unlabeled RNA per {p_key} to {p_e_path}')
        if control:
            p_e_path = estimation.estimate_p_e_control(
                df_counts,
                p_e_path,
                conversions=conversions,
            )
        elif nasc:
            rates = preprocessing.read_rates(os.path.join(count_dir, f'{constants.RATES_PREFIX}.csv'))
            if groups:
                rates['group'] = rates['barcode'].map(groups)
            p_e_path = estimation.estimate_p_e_nasc(
                rates,
                p_e_path,
                group_by=[p_key],
            )
        else:
            p_e_path = estimation.estimate_p_e(
                df_counts,
                p_e_path,
                conversions=conversions,
                group_by=[p_key],
            )
    if control:
        logger.info('Downstream processing skipped for controls')
        logger.info(f'Use `--p-e {p_e_path}` to run test samples')
        stats.end()
        stats.save(stats_path)
        return
    p_es = estimation.read_p_e(p_e_path, group_by=[p_key])

    # Aggregate counts to construct A matrix
    # NOTE: we don't use groupings here because we may need to use individual
    # barcodes later. For instance, p_c may be estimated in groups, but pi_g may
    # be estimated per cell. So that the aggregated A matrix is compatible with both
    # estimation procedures, we don't care about groupings here. Instead, groupings
    # should be manually done at each step that requires such groupings.
    aggregates_paths = {}
    for key in set(reads).union(['transcriptome'] if transcriptome_all else ['total']):
        logger.info(f'Aggregating counts for `{key}`')

        if key == 'transcriptome':
            df = df_counts[df_counts['transcriptome']]
        elif key == 'total':
            df = df_counts
        else:
            df = df_counts[df_counts['velocity'] == key]

        for convs in conversions:
            convs = sorted(convs)
            other_convs = list(set(all_conversions) - set(convs))
            aggregates_paths.setdefault(key, {})[tuple(convs)] = preprocessing.aggregate_counts(
                df[(df[other_convs] == 0).all(axis=1)] if other_convs else df,
                os.path.join(out_dir, f'A_{key}_{"_".join(convs)}.csv'),
                conversions=convs,
            )

    # Estimate p_c
    p_c_paths = {}
    for convs in conversions:
        convs = sorted(convs)
        p_c_path = os.path.join(out_dir, f'{constants.P_C_PREFIX}_{"_".join(convs)}.csv')
        logger.info(f'Estimating {convs} conversion rate in labeled RNA per {p_key} to {p_c_path}')
        df_aggregates = preprocessing.read_aggregates(
            aggregates_paths['transcriptome' if transcriptome_all else 'total'][tuple(convs)]
        )
        if groups:
            df_aggregates['group'] = df_aggregates['barcode'].map(groups).astype('category')

        p_c_paths[tuple(convs)] = estimation.estimate_p_c(
            df_aggregates, p_es, p_c_path, group_by=[p_key], threshold=cell_threshold, n_threads=n_threads, nasc=nasc
        )
    p_cs = {tuple(convs): estimation.read_p_c(p_c_paths[tuple(convs)], group_by=[p_key]) for convs in conversions}

    logger.info(f'Compling STAN model from {config.MODEL_PATH}')
    model = pystan.StanModel(file=config.MODEL_PATH, model_name=config.MODEL_NAME)

    pis = None
    alphas = None
    if method == 'pi_g':
        pi_key = 'barcode' if ignore_groups_for_est or not groups else 'group'
        pi_paths = {}
        pi_as = {}
        pi_bs = {}
        pis = {}
        for key in reads:
            for convs in conversions:
                convs = sorted(convs)
                pi_path = os.path.join(out_dir, f'pi_{key}_{"_".join(convs)}.csv')
                logger.info(
                    f'Estimating fraction of labeled `{key}` RNA for conversions {convs} per {pi_key}-gene to {pi_path}'
                )
                df_aggregates = preprocessing.read_aggregates(aggregates_paths[key][tuple(convs)])
                if groups:
                    df_aggregates['group'] = df_aggregates['barcode'].map(groups).astype('category')

                pi_paths.setdefault(key, {})[tuple(convs)] = estimation.estimate_pi(
                    df_aggregates,
                    p_es,
                    p_cs[tuple(convs)],
                    pi_path,
                    group_by=[pi_key, 'GX'],
                    p_group_by=[p_key],
                    n_threads=n_threads,
                    threshold=cell_gene_threshold,
                    seed=seed,
                    nasc=nasc,
                    model=model,
                )
                pi_a, pi_b, pi = estimation.read_pi(pi_path, group_by=[pi_key, 'GX'])
                pi_as.setdefault(key, {})[tuple(convs)] = pi_a
                pi_bs.setdefault(key, {})[tuple(convs)] = pi_b
                pis.setdefault(key, {})[tuple(convs)] = pi

        # Estimated pis need to be per cell because the adata is per cell
        if groups and not ignore_groups_for_est:
            pi_as = {
                key: {
                    convs: {(barcode, gx): value
                            for (group, gx), value in pi_as[key][convs].items()
                            for barcode in group_cells[group]}
                    for convs in pi_as[key]
                }
                for key in pi_as
            }
            pi_bs = {
                key: {
                    convs: {(barcode, gx): value
                            for (group, gx), value in pi_bs[key][convs].items()
                            for barcode in group_cells[group]}
                    for convs in pi_bs[key]
                }
                for key in pi_bs
            }
            pis = {
                key: {
                    convs: {(barcode, gx): value
                            for (group, gx), value in pis[key][convs].items()
                            for barcode in group_cells[group]}
                    for convs in pis[key]
                }
                for key in pis
            }
    elif method == 'alpha':
        pi_c_key = 'group' if groups else 'barcode'
        alpha_key = 'group' if groups and not ignore_groups_for_est else 'barcode'
        pi_c_paths = {}
        pi_cs = {}
        pi_as = {}
        pi_bs = {}
        alpha_paths = {}
        alphas = {}
        for key in reads:
            for convs in conversions:
                convs = sorted(convs)
                pi_c_path = os.path.join(out_dir, f'pi_{key}_{"_".join(convs)}.csv')
                logger.info(
                    f'Estimating fraction of labeled `{key}` RNA for conversions {convs} per {pi_c_key} to {pi_c_path}'
                )
                df_aggregates = preprocessing.read_aggregates(aggregates_paths[key][tuple(convs)])
                if groups:
                    df_aggregates['group'] = df_aggregates['barcode'].map(groups).astype('category')

                pi_c_path = estimation.estimate_pi(
                    df_aggregates,
                    p_es,
                    p_cs[tuple(convs)],
                    pi_c_path,
                    group_by=[pi_c_key],
                    p_group_by=[p_key],
                    n_threads=n_threads,
                    threshold=cell_gene_threshold,
                    seed=seed,
                    nasc=nasc,
                    model=model,
                )
                pi_c_paths.setdefault(key, {})[tuple(convs)] = pi_c_path

                pi_as.setdefault(key, {})[tuple(convs)], pi_bs.setdefault(key,
                                                                          {})[tuple(convs)], pi_c = estimation.read_pi(
                                                                              pi_c_path, group_by=[pi_c_key]
                                                                          )
                pi_cs.setdefault(key, {})[tuple(convs)] = pi_c

                alpha_path = os.path.join(out_dir, f'alpha_{key}_{"_".join(convs)}.csv')
                logger.info(
                    f'Estimating detection rate of `{key}` RNA for conversions {convs} per {alpha_key} to {alpha_path}'
                )
                alpha_path = estimation.estimate_alpha(
                    df_counts,
                    pi_c,
                    alpha_path,
                    conversions=convs,
                    group_by=[alpha_key],
                    pi_c_group_by=[pi_c_key],
                )
                alpha_paths.setdefault(key, {})[tuple(convs)] = alpha_path
                alphas.setdefault(key, {})[tuple(convs)] = estimation.read_alpha(alpha_path, group_by=[alpha_key])
        if groups and not ignore_groups_for_est:
            alphas = {
                key: {
                    convs:
                    {barcode: value
                     for group, value in alphas[key][convs].items()
                     for barcode in group_cells[group]}
                    for convs in alphas[key]
                }
                for key in alphas
            }

    else:
        raise Exception(f'Unrecognized method {method}')

    adata_path = os.path.join(out_dir, constants.ADATA_FILENAME)
    logger.info(f'Combining results into Anndata object at {adata_path}')
    adata = utils.results_to_adata(
        df_counts, conversions, gene_infos=gene_infos if not by_name else None, pis=pis, alphas=alphas
    )
    # If groups were provided, add the group as a column
    if groups:
        adata.obs['group'] = adata.obs.index.map(groups).astype('category')
    # Add the count dir that was provided as input as another column if multiple
    # count dirs were provided
    if len(count_dirs) > 1:
        adata.obs['count_dir'] = adata.obs.index.str.split('-').str[-1].astype(int).map({
            i: count_dir
            for i, count_dir in enumerate(count_dirs)
        }).astype('category')

    adata.obs.reset_index(inplace=True)

    # Add p_e, p_c estimates
    adata.obs['p_e'] = adata.obs[p_key].map(p_es).astype(float)
    for convs in conversions:
        convs = sorted(convs)
        convs_key = "_".join(convs)
        adata.obs[f'p_c_{convs_key}'] = adata.obs[p_key].map(p_cs[tuple(convs)]).astype(float)

    adata.obs.set_index('barcode', inplace=True)

    adata.write(adata_path, compression='gzip')
    stats.end()
    stats.save(stats_path)
