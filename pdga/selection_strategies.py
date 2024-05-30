import numpy as np 
from random import sample

def ranking(scores:list, k:int):
    """
    Returns the indices of the top k ranking sequences. 
    """
    idx = np.argsort(scores)[:k]
    return idx

def random(scores:list, k:int):
    """
    Returns the indices of k randomly selected sequences.
    """
    idxs = [i for i in range(0, len(scores))]
    idx = sample(idxs, k)
    return idx

def get_selection_strategy(selection_strategy:str):
    """
    Returns the function of the selected selection strategy.
    """
    if selection_strategy not in ['ranking', 'random']:
        raise ValueError('Selection strategy should be "ranking" or "random')
    else:
        if selection_strategy == 'ranking':
            return ranking
        elif selection_strategy == 'random':
            return random 