import re
from random import choice, choices, randint

# Following functions contain the basic string operations to add, remove, mutate and cyclize any sequences. 

def bb_deletion(seq:str):
    """
    Given a sequence string, removes a random occurrence of 'BB' followed by 3 digits.
    """
    matches = list(re.finditer(r'BB[0-9]{3}-', seq))
    replace = choice(matches)
    seq_new = seq[:replace.start()] + '' + seq[replace.end():]
    return seq_new

def bb_insertion(seq:str, bb_list:list):
    """
    Given a sequence string, replaces a random occurrence of '-' with a new building block.
    """
    matches = list(re.finditer(r'-', seq))
    replace = choice(matches)
    seq_new = seq[:replace.start()] + f'-{choice(bb_list)}-' + seq[replace.end():]
    return seq_new

def bb_replacement(seq:str, bb_list:list):
    """
    Given a sequence string, replaces a random occurrence of 'BB' followed by three digits 
    with a new building block.
    """
    matches = list(re.finditer(r'BB[0-9]{3}', seq))
    replace = choice(matches)
    seq_new = seq[:replace.start()] + choice(bb_list) + seq[replace.end():]
    return seq_new

def ncap_deletion(seq:str):
    """
    Given a sequence string, removes a random occurrence of 'T' followed by 3 digits.
    """
    matches = list(re.finditer(r'T[0-9]{3}-', seq))
    replace = choice(matches)
    seq_new = seq[:replace.start()] + '' + seq[replace.end():]
    return seq_new

def ncap_insertion(seq:str, ncap_list:list):
    """
    Given a sequence string, replaces a random occurrence of '-' with a new N-cap.
    """
    seq_new = f'{choice(ncap_list)}-' + seq
    return seq_new

def ncap_replacement(seq:str, ncap_list:list):
    """
    Given a sequence string, replaces a random occurrence of 'T' followed by three digits 
    with a new building block.
    """
    matches = list(re.finditer(r'T[0-9]{3}', seq))
    replace = choice(matches)
    seq_new = seq[:replace.start()] + choice(ncap_list) + seq[replace.end():]
    return seq_new

def branch_deletion(seq:str):
    """
    Given a sequence string, removes a random occurrence of 'b' followed by 3 digits.
    """
    matches = list(re.finditer(r'b[0-9]{3}-', seq))
    replace = choice(matches)
    seq_new = seq[:replace.start()] + '' + seq[replace.end():]
    return seq_new

def branch_insertion(seq:str, branch_list:list):
    """
    Given a sequence string, replaces a random occurrence of '-' with a new branching residue.
    """
    matches = list(re.finditer(r'-', seq))
    replace = choice(matches)
    seq_new = seq[:replace.start()] + f'-{choice(branch_list)}-' + seq[replace.end():]
    return seq_new

def branch_replacement(seq:str, branch_list:list):
    """
    Given a sequence string, replaces a random occurrence of 'b' followed by three digits 
    with a new branching residue.
    """
    matches = list(re.finditer(r'b[0-9]{3}', seq))
    replace = choice(matches)
    seq_new = seq[:replace.start()] + choice(branch_list) + seq[replace.end():]
    return seq_new

def c_to_n_cyclization(seq:str):
    """
    Given a sequence string, adds the cyclization term 'c' at the start and end of the sequence.
    """
    return f'c-{seq}-c'

def break_c_to_n_cyclization(seq:str):
    """
    Given a sequence string, removes the cyclization term 'c' at the start and end of the sequence.
    """
    return re.sub(r'-?c-?', '', seq)

def disulfide_full_cyclization(seq:str):
    """
    Given a sequence string, adds the cyclization term 's0' at start and end of the sequence.
    """
    return f's0-{seq}-s0'

def break_disulfide_full_cyclization(seq:str):
    """
    Given a sequence string, removes the cyclization term 's0' at the start and end of the sequence.
    """
    return re.sub(r'-?s0-?', '', seq)

def disulfide_cyclization(seq:str):
    """
    Given a sequence string, adds the cyclization term 'sx' at two random positions in the sequence.
    """
    if 's1' not in seq:
        match_1 = list(re.finditer(r'-', seq))
        replace_1 = choice(match_1)
        seq_1 = seq[:replace_1.start()] +  f'-s1-' + seq[replace_1.end():] 
        match_2 = list(re.finditer(r'-', seq_1))
        replace_2 = choice(match_2)
        seq_new = seq_1[:replace_2.start()] +  f'-s1-' + seq_1[replace_2.end():]
        return seq_new
    elif 's2' not in seq:
        match_1 = list(re.finditer(r'-', seq))
        replace_1 = choice(match_1)
        seq_1 = seq[:replace_1.start()] +  f'-s2-' + seq[replace_1.end():] 
        match_2 = list(re.finditer(r'-', seq_1))
        replace_2 = choice(match_2)
        seq_new = seq_1[:replace_2.start()] +  f'-s2-' + seq_1[replace_2.end():]
        return seq_new
    elif 's3' not in seq:
        match_1 = list(re.finditer(r'-', seq))
        replace_1 = choice(match_1)
        seq_1 = seq[:replace_1.start()] +  f'-s3-' + seq[replace_1.end():] 
        match_2 = list(re.finditer(r'-', seq_1))
        replace_2 = choice(match_2)
        seq_new = seq_1[:replace_2.start()] +  f'-s3-' + seq_1[replace_2.end():]
        return seq_new
    else:
        return seq

def break_disulfide_cyclization(seq:str):
    """
    Given a sequence string, removes the cyclization term 'sx' in the sequence.
    """
    if ('s1' in seq) and ('s2' in seq) and ('s3' in seq):
        return re.sub('s3-', '', seq)
    elif ('s1' in seq) and ('s2' in seq):
        return re.sub('s2-', '', seq)
    elif 's1' in seq:
        return re.sub('s1-', '', seq)
    else:
        return seq

# Following functions combine the previous basic functions into functional subunits. 

def bb_mutations(seq:str, bb_list:list):
    """
    Combines building block deletion, insertion and replacement functions. 
    Chooses and applies one randomly selected mutation. 
    """
    mutation_type = choice(['deletion',
                            'insertion',
                            'replacement'
                            ])

    if mutation_type == 'deletion':
        return bb_deletion(seq)
    elif mutation_type == 'insertion':
        return bb_insertion(seq, bb_list)
    elif mutation_type == 'replacement':
        return bb_replacement(seq, bb_list)
    else:
        return seq

def ncap_mutations(seq:str, ncap_list:list):
    """
    Combines N-cap deletion, insertion and replacement functions. 
    If the sequence contains an ncap, it chooses and applies deletion or replacement of the N-cap. 
    If the sequence does not contain an ncap and doesn't contain a cyclization term, it applies insertion of an N-cap.
    """
    if 'T' in seq:
        mutation_type = choice(['deletion', 'replacement'])
        if mutation_type == 'deletion':
            return ncap_deletion(seq)
        elif mutation_type == 'replacement':
            return ncap_replacement(seq, ncap_list)
    elif ('T' not in seq) and ('c' not in seq) and ('s0' not in seq):
        return ncap_insertion(seq, ncap_list)
    else:
        return seq
    
def branch_mutations(seq:str, branch_list:list):
    """
    Combines branching residue deletion, insertion and replacement functions. 
    If the sequence contains a branching residue, it chooses and applies deletion or replacement of the branch. 
    If the sequence does not contain a branching residue and doesn't contain a cyclization term, it applies insertion of a branch.
    """
    if 'b' in seq:
        mutation_type = choice(['deletion', 'replacement'])
        if mutation_type == 'deletion':
            return branch_deletion(seq)
        elif mutation_type == 'replacement':
            return branch_replacement(seq, branch_list)
    elif ('b' not in seq) and ('c' not in seq) and ('s0' not in seq):
        return branch_insertion(seq, branch_list)
    else:
        return seq

def cyclization_mutations(seq:str):
    """
    Combines all cyclization functions. 
    If the sequence does not contain N-caps, branching residues, c-to-n-cyclizations or disulfide cyclization, a c-to-n-cyclization is inserted.
    If the sequence contains a c-to-n-cyclization term, it removes the cyclization. 
    """
    if ('c' not in seq) and ('b' not in seq) and ('T' not in seq) and ('s0' not in seq):
        return c_to_n_cyclization(seq)
    elif 'c' in seq:
        return break_c_to_n_cyclization(seq)
    else:
        return seq
    
def disulfide_mutations(seq:str):
    """
    Combines disulfide cyclization and disulfide full cyclization functions. 
    If the sequence does not contain a disulfide cyclization term, it inserts a disulfide cyclization. 
    If the sequence does not contain a disulfide full cyclization term, it inserts a disulfide full cyclization. 
    If the sequence contains a disulfide cyclization term, it removes the cyclization. 
    If the sequence contains a disulfide full cyclization term, it removes the cyclization.
    """
    if ('s' not in seq) and ('b' not in seq) and ('T' not in seq) and ('c' not in seq):
        mutation_type = choice(['disulfide_cyclization', 'disulfide_full_cyclization'])
        if mutation_type == 'disulfide_cyclization':
            return disulfide_cyclization(seq)
        elif mutation_type == 'disulfide_full_cyclization':
            return disulfide_full_cyclization(seq)
    elif ('s0' not in seq) and ('s' in seq) and ('b' not in seq) and ('T' not in seq) and ('c' not in seq):
        mutation_type = choice(['disulfide_full_cyclization', 'break_disulfide_cyclization'])
        if mutation_type == 'disulfide_full_cyclization':
            return disulfide_full_cyclization(seq)
        elif mutation_type == 'break_disulfide_cyclization':
            return break_disulfide_cyclization(seq)
    elif ('s0' in seq) and ('s' in seq):
        mutation_type = choice(['break_disulfide_full_cyclization', 'break_disulfide_cyclization'])
        if mutation_type == 'break_disulfide_full_cyclization':
            return break_disulfide_full_cyclization(seq)
        elif mutation_type == 'break_disulfide_cyclization':
            return break_disulfide_cyclization(seq)
    elif ('s0' not in seq) and ('s' in seq):
        return break_disulfide_cyclization(seq)
    elif 's0' in seq:
        return break_disulfide_full_cyclization(seq)
    else:
        return seq

def fixed_positions_mutations(seq:str, bb_list:list, ncap_list:list, branch_list:list):
    """
    Combines building block, N-cap and branching residue replacement. 
    Checks which of these are found in the provided template and replaces them without changing the building block type or 
    the length of the template sequence to insert. 

    """
    bb_types = []
    if 'BB' in seq:
        bb_types.append('bb')
    if 'T' in seq:
        bb_types.append('ncap')
    if 'b' in seq:
        bb_types.append('branch')

    mutation_type = choice(bb_types)
    if mutation_type == 'bb':
        return bb_replacement(seq, bb_list)
    elif mutation_type == 'ncap':
        return ncap_replacement(seq, ncap_list)
    elif mutation_type == 'branch':
        return branch_replacement(seq, branch_list)
    else:
        return seq

#Following functions provide the interaction between the GA and the mutation functions defined before.

def topology_to_mutation_types(topology:str):
    """
    Returns the type of allowed mutations and their cumulative weights based on the selected topology. 
    Changing the weights will change the frequency with which a type of mutation will occur. 
    """
    if topology == 'free':
        mutations = ['building block', 'ncap', 'branch', 'cyclization']
        weights = [70, 10, 10, 10]
    elif topology == 'free_disulfide':
        mutations = ['building block', 'ncap', 'branch', 'cyclization', 'disulfide']
        weights = [70, 5, 5, 10, 10]
    elif topology == 'linear':
        mutations = ['building block', 'ncap']
        weights = [90, 10]
    elif topology == 'cyclic':
        mutations = ['building block', 'cyclization']
        weights = [75, 25]
    elif topology == 'fixed positions':
        mutations = ['fixed positions']
        weights = [100]
    return mutations, weights

def mutate(seq:str, bb_list:list, ncap_list:list, branch_list:list, mutations:list, weights:list):
    """
    Applies one random mutation based on the provided mutations list and their cumulative weights.
    """
    mutation_type = choices(population = mutations, weights = weights, k = 1)[0]
    seq_list = seq.split('-')
    if len(seq_list) > 3:
        if mutation_type == 'building block':
            return bb_mutations(seq, bb_list)
        elif mutation_type == 'ncap':
            return ncap_mutations(seq, ncap_list)
        elif mutation_type == 'branch':
            return branch_mutations(seq, branch_list)
        elif mutation_type == 'cyclization':
            return cyclization_mutations(seq)
        elif mutation_type == 'disulfide':
            return disulfide_mutations(seq)
        elif mutation_type == 'fixed positions':
            return fixed_positions_mutations(seq, bb_list, ncap_list, branch_list)
    else:
        return seq

def crossover(parents:tuple):
    """
    Splits two parent sequences into two sub sequences and recombines them cross-wise.
    """
    parent1, parent2 = parents
    if (('c' in parent1) and ('s' in parent1)) or (('c' in parent2) and ('s' in parent2)):
        match1 = list(re.finditer(r'-', parent1))
        match2 = list(re.finditer(r'-', parent2))
        replace1 = choice(match1)
        replace2 = choice(match2)
        seq_new1 = parent1[:replace1.start()] + f'-' + parent2[replace2.end():]
        seq_new2 = parent2[:replace2.start()] + f'-' + parent1[replace1.end():]
        return seq_new1, seq_new2
    else:
        return parent1, parent2