import os
from datetime import date 
import random
from random import choice, sample
import numpy as np
import pandas as pd
from rdkit import Chem
from tqdm import tqdm

from .sequence_operations import process_building_blocks, random_linear_seq, template_to_smiles, seq_to_smiles, count_aa_and_peptoids_in_mol
from .selection_strategies import get_selection_strategy
from .fingerprints import get_fingerprint
from .fitness_functions import get_fitness_function
from .mutations import topology_to_mutation_types, mutate, crossover

class PDGA:
  def __init__(self, 
               query, 
               query_format='smiles', 
               topology='free',
               template=None, 
               pop_size=50, 
               pop_selection=10, 
               mutation_ratio=0.75, 
               cut_off=0.5, 
               selection_strategy='ranking', 
               descriptor='MAP4', 
               n_iterations=1000, 
               run_id='pdga_run', 
               seed=0,
               verbose=True
               ):

    # Verbosity 
    self.verbose = not verbose

    # Set random seed.
    self.seed = seed
    random.seed(self.seed)
    np.random.seed(self.seed)

    # Determine desired descriptor. 
    self.descriptor = descriptor
    self.fingerprint = get_fingerprint(self.descriptor)
    self.fitness_func = get_fitness_function(self.descriptor)

    # Import building block lists and translation dictionary. Check if lists are compatible with the selected topology.
    bb_list, ncap_list, branch_list, translation_dict = process_building_blocks()

    if len(bb_list) == 0:
      raise ValueError('Empty building block list')
    else:
      self.bb_list = bb_list
    if len(ncap_list) == 0:
      self.ncap_list = False
    else:
      self.ncap_list = ncap_list
    if len(branch_list) == 0:
      self.branch_list = False
    else:
      self.branch_list = branch_list

    self.translation_dict = translation_dict

    # Check query format and translate query to selected fingerprint.
    if query_format == 'smiles' and descriptor == 'MAP4_COMBO':
      query_smiles = query.split('.')
      self.query_mol = [Chem.MolFromSmiles(smiles) for smiles in query_smiles]
      self.query = [self.fingerprint(mol) for mol in self.query_mol]
    elif query_format == 'smiles':
      self.query_mol = Chem.MolFromSmiles(query)
      self.query = self.fingerprint(self.query_mol)
    elif query_format == 'sequence':
      self.query_mol = Chem.MolFromSequence(query)
      self.query = self.fingerprint(self.query_mol)
    elif query_format == 'mol':
      self.query_mol = query
      self.query = self.fingerprint(self.query_mol)
    elif query_format == 'advanced_sequence':
      self.query_mol = Chem.MolFromSmiles(seq_to_smiles(query, self.translation_dict))
      self.query = self.fingerprint(self.query_mol)
    elif query_format == 'fingerprint':
      self.query = query
    else:
      raise ValueError('Query format should be one of the following: "smiles", "mol", "sequence", "advanced sequence" or "fingerprint"')

    # Determine topology used for the structure generation.
    if topology not in ['free', 'linear', 'cyclic', 'fixed positions']:
      raise ValueError('Topology should be one of the following: "free", "linear", "cyclic" or "fixed positions"')
    else:
      self.topology = topology

    # If the PDGA is used with fixed positions, extract template scaffold and template fill. 
    if template:
      self.template_scaffold = template[0]
      self.template_fill = template[1]
    
    # Get selected selection strategy.
    self.selection_strategy = get_selection_strategy(selection_strategy)

    # Initiate remaining variables and parameters.
    self.pop_size = pop_size
    self.pop_selection = pop_selection
    self.mutation_ratio = mutation_ratio
    self.cut_off = cut_off
    self.n_iterations = n_iterations
    self.run_id = run_id
    self.date = date.today().strftime('%Y%m%d')

    # Determine types of allowed mutations according to selected topology. Remove N-caps and branching mutations if building block lists are empty. 
    self.mutations, self.weights = topology_to_mutation_types(self.topology)

    if (self.ncap_list == False) and ('ncap' in self.mutations):
      idx = self.mutations.index('ncap')
      self.mutations.pop(idx)
      self.weights.pop(idx)
    
    if (self.branch_list == False) and ('branch' in self.mutations):
      idx = self.mutations.index('branch')
      self.mutations.pop(idx)
      self.weights.pop(idx)
    
    # Create starting population. With fixed positions, the starting population will be the template fill.
    if self.topology == 'fixed positions':
      self.pop_sequences = [self.template_fill for i in range(self.pop_size)]
      self.pop_smiles = [template_to_smiles(self.template_scaffold, seq, self.translation_dict) for seq in self.pop_sequences]
      self.pop_mols = [Chem.MolFromSmiles(smiles) for smiles in self.pop_smiles]
      self.scores = [self.fitness_func(self.query, mol) for mol in self.pop_mols]
    else:
      self.pop_sequences = [random_linear_seq(self.bb_list, min_len=5, max_len=30) for i in range(self.pop_size)]
      self.pop_smiles = [seq_to_smiles(seq, self.translation_dict) for seq in self.pop_sequences]
      self.pop_mols = [Chem.MolFromSmiles(smiles) for smiles in self.pop_smiles]
      self.scores = [self.fitness_func(self.query, mol) for mol in self.pop_mols]

    # Generate output data frame and initiate generation.  
    self.output = pd.DataFrame(zip(self.pop_sequences, self.pop_smiles, self.scores), columns=['sequence', 'smiles', 'dist'])
    self.generations = pd.DataFrame(columns=['generation', 'best_score', 'best_smiles', 'best_sequence'])
  
    # Store parameters in log file.
    self.parameters = locals()

    if not os.path.exists('results'):
      os.makedirs('results')

    if not os.path.exists('results/logs'):
      os.makedirs('results/logs')

    with open(f'results/logs/{self.date}_{self.run_id}_{self.descriptor}_{self.seed}_log.txt', 'w') as f:
      for key, value in self.parameters.items():
        f.write('%s:%s\n' % (key, value))
      
    if not os.path.exists('results/generations'):
      os.makedirs('results/generations')

    self.early_stop = False
    self.early_stop_counter = 0

 # Main function of the genetic algorithm. 
  def optimize(self):

    for i in tqdm(range(self.n_iterations), disable=self.verbose):

      # Select indices of generated sequences according to defined strategy. 
      selected_idx = self.selection_strategy(self.scores, self.pop_selection)
      
      # Retrieve selected sequences with the previously obtained indices. 
      selected_population = [self.pop_sequences[idx] for idx in selected_idx]
      
      # Mutations and crossovers. 
      mutated_population = []
      mut_iter = 0
      while mut_iter < int(self.mutation_ratio*(self.pop_size - self.pop_selection)):
        new_seq = mutate(choice(selected_population), self.bb_list, self.ncap_list, self.branch_list, self.mutations, self.weights)
        mutated_population.append(new_seq)
        mut_iter += 1

      crossover_population = []
      cross_iter = 0
      while mut_iter < int((1-self.mutation_ratio)*(self.pop_size - self.pop_selection)):
        new_seq1, new_seq2 = crossover(sample(selected_population, 2))
        crossover_population.extend((new_seq1, new_seq2))
        cross_iter += 1

      # Add new sequences to the working population.
      if self.topology == 'fixed positions':
        self.pop_sequences = selected_population + mutated_population + crossover_population
        self.pop_smiles = [template_to_smiles(self.template_scaffold, seq, self.translation_dict) for seq in self.pop_sequences]
        self.pop_mols = [Chem.MolFromSmiles(smiles) for smiles in self.pop_smiles]
        self.scores = [self.fitness_func(self.query, mol) for mol in self.pop_mols]
      else:
        self.pop_sequences = selected_population + mutated_population + crossover_population
        self.pop_smiles = [seq_to_smiles(seq, self.translation_dict) for seq in self.pop_sequences]
        self.pop_mols = [Chem.MolFromSmiles(smiles) for smiles in self.pop_smiles]
        self.scores = [self.fitness_func(self.query, mol) for mol in self.pop_mols]

      # Early stopping if the best score is 0 for 5 consecutive iterations.
      if any(score == 0 for score in self.scores):
        self.early_stop = True
      
      if self.early_stop:
        self.early_stop_counter += 1

      if self.early_stop_counter == 5:
        print('Early stopping due to 5 consecutive iterations with a score of 0.')
        break

      # At the end of the iteration, remove sequences with distances lower than the defined cut-off value and append to output data frame.
      df_iter = pd.DataFrame(zip(self.pop_sequences, self.pop_smiles, self.scores), columns=['sequence', 'smiles', 'dist'])
      df_iter = df_iter[df_iter['dist'] <= self.cut_off].copy()
      df_iter['generation'] = i
      self.output = pd.concat((self.output, df_iter))

      # To the generations data frame, append the best sequence, smiles and score for the current iteration.
      best_idx = self.scores.index(min(self.scores))
      self.generations = self.generations.append({'generation': i, 'best_score': self.scores[best_idx], 'best_smiles': self.pop_smiles[best_idx], 'best_sequence': self.pop_sequences[best_idx]}, ignore_index=True)
      self.generations.to_csv(f'results/generations/{self.date}_{self.run_id}_{self.descriptor}_{self.seed}.csv', index=False)
 
    # When n_iterations are finished, remove duplicates and generate the ouput csv file containing sequence, SMILES and distance to query.
    self.output = self.output.sort_values(by=['dist', 'generation'], ignore_index=True)
    self.output = self.output.drop_duplicates(subset=['sequence'], keep='first')
    self.output.to_csv(f'results/{self.date}_{self.run_id}_{self.descriptor}_{self.seed}.csv', index=False)

    