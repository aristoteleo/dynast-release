.PHONY : install test check build docs clean push_release

test:
	rm -f .coverage
	pytest -vv tests/ --cov=dynast
	coverage combine --append

check:
	flake8 dynast && echo OK
	yapf -r --diff dynast && echo OK

build:
	python setup.py sdist

docs:
	sphinx-build -a docs docs/_build

clean:
	rm -rf build
	rm -rf dist
	rm -rf dynast_release.egg-info
	rm -rf docs/_build
	rm -rf docs/api

bump_patch:
	bumpversion patch

bump_minor:
	bumpversion minor

bump_major:
	bumpversion major

push_release:
	git push && git push --tags
