import configparser
from functools import partial

from constraints import Constraints, constraint_included, constraint_distances, constraint_max_charge, constraint_max_num_changes
from evolution import *
from logger import FileLogger
from utils import *

# PARSING CONFIG
config = configparser.ConfigParser()
config.read('config.ini')

pdb_file = config['PDB']['File']
value = float(config['PDB']['VALUE'])
cros_prob = float(config['PARAMS']['CrosProb'])
mut_prob = float(config['PARAMS']['MutProb'])
mut_num = int(config['PARAMS']['MutNum'])
eval_param = float(config['PARAMS']['EvalParam'])
pop_size = int(config['PARAMS']['PopSize'])
compute_lmb_inf = config['COMPUTING']['ComputeLambdaInf']
compute_lmb_ouf = config['COMPUTING']['ComputeLambdaOuf']
computed_proteins_path = config['COMPUTING']['ComputedProteinsFileName']
result_file_name = config['COMPUTING']['ResultFileName']

logger = FileLogger("logout")

# GENERATING CONSTRAINTS
constraints = Constraints()

coordinates = read_coordinates(pdb_file)
sequence = read_sequence(pdb_file)

f1 = partial(constraint_included, aminoacids_set="DE", positions_set=PositionsSet1)
f2 = partial(constraint_distances, min_distance=5.0, coords=coordinates, positions_set=PositionsSetUnion)
f3 = partial(constraint_max_charge, max_charge=7)
f4 = partial(constraint_max_num_changes, max_num_changes=10)

constraints.add(f1)
constraints.add(f2)
constraints.add(f3)
constraints.add(f4)

# COMPUTING
population = ProteinEvolution(population=None, mut_prob=mut_prob, mut_num=mut_num, cros_prob=cros_prob,
                              input_file=compute_lmb_inf, output_file=compute_lmb_ouf, save_file=computed_proteins_path,
                              logger=logger, checker=constraints)
population.load_computed_proteins()
population.generate_population(default_sequence=sequence, default_value=value, pop_size=pop_size, from_computed=True)

iteration, step, stop_step = 1, 0, 5000

the_best_value = 0
while step < stop_step:
    logger(f"Iteration: {iteration}\n")

    population.mutation(attempts=4000)
    population.crossover(attempts=4000)
    population.compute()

#    population.selection(eval_param=0.05, save_n_best=3)
    population.selection(eval_param=eval_param, save_n_best=3)

    cur_best_value = population.get_best_protein().value
    if the_best_value < cur_best_value:
        the_best_value = cur_best_value
        step = 0
    else:
        step += 1

    logger(f"Current population:\n")
    population.print_current_population()
    logger(f"The best value: {the_best_value}\n"
           f"Step/Stop {step}/{stop_step}\n")
    logger("\n")

    iteration += 1
