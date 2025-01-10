import pickle


def read_pickle(path):
    with open(path, 'rb') as fp:
        proteins = pickle.load(fp)
    return proteins
