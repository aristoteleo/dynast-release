import pickle

import pystan


def prepare_model(model_path, model_name, pickle_path):
    """
    """
    model = pystan.StanModel(file=model_path, model_name=model_name)
    with open(pickle_path, 'wb') as f:
        pickle.dump(model, f)
    return pickle_path


def load_model(pickle_path):
    """
    """
    with open(pickle_path, 'rb') as f:
        return pickle.load(f)
