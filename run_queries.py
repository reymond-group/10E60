from pdga.pdga import PDGA
from multiprocessing import Pool

# Import queries
file_path = 'queries.csv'
with open(file_path, 'r') as file:
    lines = file.readlines()
    lines = [line.strip() for line in lines if not line.startswith('#')]
queries = [tuple(line.split(',')) for line in lines]
queries = [(query[0], query[1], i) for query in queries for i in range(3)] # Create triplicates with different seeds

# Define genetic algorithm function 
def run_ga(args):
    name, query, seed = args
    ga = PDGA(
        query=query,
        query_format='smiles',
        topology='free',
        template=None,
        pop_size=50,
        pop_selection=15,
        mutation_ratio=0.5,
        selection_strategy='ranking',
        descriptor='MXFP',
        cut_off=300,
        n_iterations=10000,
        run_id=name,
        seed=seed,
        verbose=False
    )
    ga.optimize()

# Run multiprocessing loop 
with Pool(processes=12) as pool:
    pool.map(run_ga, queries)