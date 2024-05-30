from pdga.pdga import PDGA, seq_to_smiles
from rdkit import Chem
from multiprocessing import Pool

# Import queries
file_path = 'queries_seq_to_seq.csv'
with open(file_path, 'r') as file:
    lines = file.readlines()
    lines = [line.strip() for line in lines if not line.startswith('#')]
queries = [tuple(line.split(',')) for line in lines]
queries = [(query[0], query[1], query[2], i) for query in queries for i in range(3)] # Create triplicates with different seeds

# Define genetic algorithm function 
def run_ga(args):
    
    name, start, query, seed = args
    
    ga = PDGA(
        query=query,
        query_format='advanced_sequence',
        topology='free',
        template=None,
        pop_size=50,
        pop_selection=15,
        mutation_ratio=0.5,
        selection_strategy='ranking',
        descriptor='MAP4',
        cut_off=1,
        n_iterations=10000,
        run_id=name,
        seed=seed,
        verbose=False
    )

    ga.pop_sequences = [start for i in range(ga.pop_size)]
    ga.pop_smiles = [seq_to_smiles(seq, ga.translation_dict) for seq in ga.pop_sequences]
    ga.pop_mols = [Chem.MolFromSmiles(smiles) for smiles in ga.pop_smiles]
    ga.scores = [ga.fitness_func(ga.query, mol) for mol in ga.pop_mols]

    ga.optimize()

# Run multiprocessing loop 
with Pool(processes=12) as pool:
    pool.map(run_ga, queries)