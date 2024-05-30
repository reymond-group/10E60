import numpy as np
from rdkit.Chem import AllChem, DataStructs, rdchem, rdMolDescriptors

from rdkit.Chem.AtomPairs import Pairs
from mapchiral import mapchiral
from mxfp import mxfp 

mxfp2D = mxfp.MXFPCalculator()

def calculate_ecfp(mol:rdchem.Mol):
    """
    Calculates morgan fingerprints of an rdkit mol object.
    """
    return AllChem.GetMorganFingerprintAsBitVect(mol, 2)

def calculate_ap(mol:rdchem.Mol):
    """
    Calculates atom-pair fingerprints of an rdkit mol object.
    """
    ap = rdMolDescriptors.GetHashedAtomPairFingerprintAsBitVect(mol)
    fp = np.zeros((0,), dtype=np.int8)
    DataStructs.ConvertToNumpyArray(ap, fp)
    return fp.astype(np.bool)

def calculate_mqn(mol:rdchem.Mol):
    """
    Calculates molecular quantum numbers of an rdkit mol object.
    """
    return np.array(rdMolDescriptors.MQNs_(mol))

def calculate_mxfp(mol:rdchem.Mol):
    """
    Calculates macromolecule extended fingerprint of an rdkit mol object.
    """
    return mxfp2D.mxfp_from_mol(mol)

def calculate_map4(mol:rdchem.Mol):
    """
    Calculates minhashed atom-pair fingerprint of an rdkit mol object.
    """
    return mapchiral.encode(mol)

def get_fingerprint(fingerprint:str):
    """
    Returns the selected fingerprint.
    """
    if fingerprint not in ['ECFP', 'AP', 'MQN', 'MXFP', 'MAP4', 'MAP4_COMBO']:
        raise ValueError('Fitness function should be "ECFP", "AP", "MQN", "MXFP" or "MAP4"')
    else:
        if fingerprint == 'ECFP':
            return calculate_ecfp
        elif fingerprint == 'AP':
            return calculate_ap
        elif fingerprint == 'MQN':
            return calculate_mqn
        elif fingerprint == 'MXFP':
            return calculate_mxfp
        elif fingerprint == 'MAP4':
            return calculate_map4
        elif fingerprint == 'MAP4_COMBO':
            return calculate_map4