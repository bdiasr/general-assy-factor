import json
import matplotlib.pyplot as plt
import pathlib
import pickle
import time




def open_file(path='../pickles/', file_name='', formatting='.pickle', operation='rb'):
    """ Return a file with specific format either for
      reading (operation='rb') or writing (operation='wb')
    """
    return  open(str(pathlib.Path(path + file_name + formatting).absolute()), operation)