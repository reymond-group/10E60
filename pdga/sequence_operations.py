import pandas as pd
from random import randint, choice

from rdkit import Chem
from rdkit.Chem import rdchem

def process_building_blocks():
    """
    Process csv files containing building block lists:
    bb: amino acid or peptoid building blocks to be used. 
    ncaps: capping elements for the N-terminus 
    branches: possible branching building blocks (contain two amino moieties)
    additional: only used when using fixed positions and the query contains building blocks that GA is not using as possible mutations. 
    """
    building_blocks = pd.read_csv('pdga/BB/bb.csv', comment=';')
    ncaps = pd.read_csv('pdga/BB/ncaps.csv', comment=';')
    branches = pd.read_csv('pdga/BB/branches.csv', comment=';')
    additional = pd.read_csv('pdga/BB/additional.csv', comment=';')

    bb_dict = dict(zip(building_blocks.ID, building_blocks.SMILES))
    ncap_dict = dict(zip(ncaps.ID, ncaps.SMILES))
    branch_dict = dict(zip(branches.ID, branches.SMILES))
    additional_dict = dict(zip(additional.ID, additional.SMILES))

    bb_list = building_blocks.ID.values.tolist()
    ncap_list = ncaps.ID.values.tolist()
    branch_list = branches.ID.values.tolist()

    translation_dict = {**bb_dict,**ncap_dict, **branch_dict, **additional_dict}

    translation_dict['c'] = '9'
    translation_dict['s0'] = 'NC(CS4)C(=O)'
    translation_dict['s1'] = 'NC(CS5)C(=O)'
    translation_dict['s2'] = 'NC(CS6)C(=O)'
    translation_dict['s3'] = 'NC(CS7)C(=O)'

    return bb_list, ncap_list, branch_list, translation_dict

def random_linear_seq(bb_list:list, min_len:int=5, max_len:int=30):
    """
    Generate a random linear sequence from a provided building block list.
    """
    seq = choice(bb_list)
    for i in range(randint(min_len, max_len)):
        seq = seq + '-' + choice(bb_list)
    return seq

def seq_to_smiles(seq:str, translation_dict:dict):
    """
    Translates sequences to SMILES using a translation dictionary.
    """
    smiles = ''
    try:
        for element in seq.split('-'):
            smiles += translation_dict.get(element)
        if 'b' in seq:
            return smiles + '8'
        elif 'c' in seq:
            return smiles[1] + smiles[0] + smiles[2:]
        else:
            return smiles + 'O'
    except:
        print(f'Some error occurred with sequence: {seq}')
        return ''

def template_to_seq(template_scaffold:str, seq:str):
    """
    Recombines a provided template scaffold with a generated filling sequence.
    """
    seq_split = seq.split('-')
    return template_scaffold.format(*seq_split)

def template_to_smiles(template_scaffold:str, seq:str, translation_dict:dict):
    """
    Translates a template scaffold and a filling sequence to SMILES.
    """
    seq_new = template_to_seq(template_scaffold, seq)
    smiles = ''
    try:
        for element in seq_new.split('-'):
            smiles += translation_dict.get(element)
        if 'b' in seq_new:
            return smiles + '8'
        elif 'c' in seq_new:
            return smiles[1] + smiles[0] + smiles[2:]
        else:
            return smiles + 'O'
    except:
        print(f'Some error occurred with sequence: {seq_new}')
        return ''

def count_peptoids_in_mol(mol:rdchem.Mol=None):
    """
    Counts the number of peptoids in a mol object.
    """
    smarts = '[C][CH0][NH0][CH2;$(CC=O)]'
    matches = mol.GetSubstructMatches(Chem.MolFromSmarts(smarts))
    return len(matches)

def count_aa_and_peptoids_in_mol(mol):
    """
    Counts the number of amino acids and peptoids in a mol object.
    """
    bb_smarts = '[CH0](=O)[NH0,NH1]'
    c_to_n_cyclization_smarts = '[RCH0](=O)[RNH0,RNH1]'

    bb_matches = mol.GetSubstructMatches(Chem.MolFromSmarts(bb_smarts))
    if mol.HasSubstructMatch(Chem.MolFromSmarts(c_to_n_cyclization_smarts)):
        return len(bb_matches)
    else:
        return len(bb_matches) + 1 