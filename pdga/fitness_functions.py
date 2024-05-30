import numpy as np
from rdkit.Chem import DataStructs, rdchem, rdMolDescriptors
from mapchiral.mapchiral import jaccard_similarity

from .sequence_operations import count_peptoids_in_mol, count_aa_and_peptoids_in_mol
from .fingerprints import calculate_ecfp, calculate_ap, calculate_mqn, calculate_mxfp, calculate_map4

def morgan_dist(query, mol:rdchem.Mol, force_peptoids:int=0):
    """
    Calculates the Tanimoto distance between two Morgan fingerprints.
    """
    fp1 = query
    fp2 = calculate_ecfp(mol)
    distance = 1 - DataStructs.FingerprintSimilarity(fp1, fp2)
    if force_peptoids != 0:
            return distance + distance*(abs(force_peptoids-count_peptoids_in_mol(mol)))
    else:
        return distance

def ap_dist(query, mol:rdchem.Mol, force_peptoids:int=0):
    """
    Calculates the Tanimoto distance between two atom-pair fingerprints.
    """
    fp1 = query
    fp2 = calculate_ap(mol)

    intersection = np.logical_and(fp1, fp2)
    union = np.logical_or(fp1, fp2)

    distance = 1 - intersection.sum() / float(union.sum())
    if force_peptoids != 0:
            return distance + distance*(abs(force_peptoids-count_peptoids_in_mol(mol)))
    else:
        return distance

def mqn_dist(query, mol:rdchem.Mol, force_peptoids:int=0):
    """
    Calculates the Manhattan distance between two MQN fingerprints.
    """
    fp1 = query
    fp2 = calculate_mqn(mol)
    distance = sum(abs(fp1 - fp2))
    if force_peptoids != 0:
            return distance + distance*(abs(force_peptoids-count_peptoids_in_mol(mol)))
    else:
        return distance

def mxfp_dist(query, mol:rdchem.Mol):
    """
    Calculates the Manhattan distance between two macromolecule-extended fingerprints.
    """
    fp1 = query
    fp2 = calculate_mxfp(mol)
    return sum(abs(fp1 - fp2))

def map4_dist(query, mol:rdchem.Mol, force_peptoids:int=0):
    """
    Calculates the Jaccard distance between two minhashed atom-pair fingerprints.
    """
    fp1 = query
    fp2 = calculate_map4(mol)
    distance = 1.0 - jaccard_similarity(fp1, fp2)
    if force_peptoids != 0:
            return distance + distance*(abs(force_peptoids-count_peptoids_in_mol(mol)))*0.1
    else:
        return distance

def map4_combo_dist(query, mol:rdchem.Mol):
    """
    Calculates the Jaccard distance between two minhashed atom-pair fingerprints. Adds penalty 
    """
    fp1 = query[0]
    fp2 = query[1]
    fp3 = calculate_map4(mol)
    distance1 = 1.0 - jaccard_similarity(fp1, fp3)
    distance2 = 1.0 - jaccard_similarity(fp2, fp3)
    return distance1 + distance2 + abs(distance1 - distance2)


def get_fitness_function(fitness_function:str):
    """
    Returns the selected fitness function.
    """
    if fitness_function not in ['ECFP', 'AP', 'MQN', 'MXFP', 'MAP4', 'MAP4_COMBO']:
        raise ValueError('Fitness function should be "ECFP", "AP", "MQN", "MXFP", "MAP4" or "MAP4_COMBO"')
    else:
        if fitness_function == 'ECFP':
            return morgan_dist
        elif fitness_function == 'AP':
            return ap_dist
        elif fitness_function == 'MQN':
            return mqn_dist
        elif fitness_function == 'MXFP':
            return mxfp_dist
        elif fitness_function == 'MAP4':
            return map4_dist
        elif fitness_function == 'MAP4_COMBO':
            return map4_combo_dist