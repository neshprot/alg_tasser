import configparser
from functools import partial

from constraints import Constraints, constraint_distances, constraint_max_charge, constraint_max_num_changes
from evolution import *
from logger import FileLogger
from utils import *

# PARSING CONFIG
config = configparser.ConfigParser()
config.read('config.ini')

# задаём некоторые константы из config
pdb_file = config['PDB']['File']
value = float(config['PDB']['VALUE'])
pka215 = float(config['PDB']['pka215'])
pka207 = float(config['PDB']['pka207'])
cros_prob = float(config['PARAMS']['CrosProb'])
mut_prob = float(config['PARAMS']['MutProb'])
mut_num = int(config['PARAMS']['MutNum'])
eval_param = float(config['PARAMS']['EvalParam'])
pop_size = int(config['PARAMS']['PopSize'])
weightval = float(config['PARAMS']['Weightval'])
weight215_207 = float(config['PARAMS']['Weight215_207'])
compute_lmb_inf = config['COMPUTING']['ComputeLambdaInf']
compute_lmb_ouf = config['COMPUTING']['ComputeLambdaOuf']
computed_proteins_path = config['COMPUTING']['ComputedProteinsFileName']
result_file_name = config['COMPUTING']['ResultFileName']

logger = FileLogger("logout")

# GENERATING CONSTRAINTS
constraints = Constraints()

coordinates = read_coordinates(pdb_file)
sequence = read_sequence(pdb_file)

# функции ограничений
f1 = partial(constraint_distances, min_distance=5.0, coords=coordinates, positions_set=PositionsSetUnion)
f2 = partial(constraint_max_charge, max_charge=7)
f3 = partial(constraint_max_num_changes, max_num_changes=6)

constraints.add(f1)
constraints.add(f2)
constraints.add(f3)

# COMPUTING
population = ProteinEvolution(population=None, mut_prob=mut_prob, mut_num=mut_num, cros_prob=cros_prob,
                              input_file=compute_lmb_inf, output_file=compute_lmb_ouf, save_file=computed_proteins_path,
                              logger=logger, checker=constraints, weightval=weightval,
                              weight215_207=weight215_207)
population.load_computed_proteins()
population.generate_population(default_sequence=sequence, default_value=value, default_pka215=pka215,
                               default_pka207=pka207, pop_size=pop_size, from_computed=True)
population.generate_start_compute(default_value=value)
