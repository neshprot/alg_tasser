import math
import os
import random
import time
from abc import abstractmethod, ABC
from copy import copy

from data import Protein, Gene, POLAR, CHARGED, BULKY

Pull = "ARNDVHGQEILKMPSYTWFV"   # список 20 существующих аминокислот

PositionsSet1 = {23, 30, 49, 52, 53, 56, 60, 88, 89, 92, 93, 96, 188, 214, 215, 216, 217, 218, 220, 221, 222, 223, 224}  # within 5 of resid 219
PositionsSet2 = {23, 56, 60, 63, 85, 89, 184, 188, 210, 211, 212, 213, 214, 216, 217, 218, 219, 220, 221} # within 5 of resid 215
PositionsSet3 = {1, 80, 85, 191, 192, 194, 195, 196, 197, 198, 202, 203, 204, 205, 206, 208, 209, 210, 211, 212}
PositionsSetUnion = set.union(PositionsSet1, PositionsSet2, PositionsSet3)

class Evolution(ABC):
    def __init__(self):
        self._population = None
        self._mut_prob = None
        self._cros_prob = None

    @property
    def population(self):
        return self._population

    @property
    def mut_prob(self):
        return self._mut_prob

    @property
    def cros_prob(self):
        return self._cros_prob

    @abstractmethod
    def mutation(self, *args, **kwargs):
        pass

    @abstractmethod
    def crossover(self, *args, **kwargs):
        pass

    @abstractmethod
    def selection(self, *args, **kwargs):
        pass


class BaseFunction(ABC):
    def __init__(self):
        self._input_file: str = None
        self._output_file: str = None
        self._save_file: str = None
        self._computed = None

    @abstractmethod
    def compute(self, *args, **kwargs):
        pass

    @abstractmethod
    def save_computing(self, *args, **kwargs):
        pass

# класс с основными функциями эфолюции
class ProteinEvolution(Evolution, BaseFunction):
    def __init__(self, population, mut_prob, mut_num, cros_prob, input_file,
                 output_file, save_file, logger, checker=None,
                 weightval=1.0, weight215_207=1.0):
        super().__init__()
        self._population = population
        self._mut_prob = mut_prob
        self._mut_num = mut_num
        self._cros_prob = cros_prob
        self._input_file = input_file
        self._output_file = output_file
        self._save_file = save_file
        self._computed = dict()
        self._checker = checker
        self._logger = logger
        self._weightval = weightval
        self._weight215_207 = weight215_207

    # nutation function
    def mutation(self, attempts=1, iteration=1, sets={}, pulls={}, probs={}):
        """

        :param attempts: число попыток инциниализации на один protein
        :return:
        """

        new_population = []     # список белков в новой популяции
        num_of_changed = 0      # кол-во измененных белков

        # перебор белков в популяции
        for protein in self._population:
            new_protein = copy(protein)
            # условие возникновения мутации(с вероятностью mut_prob)
            if random.random() < self._mut_prob:
                attempt = 0
                num_of_changed += 1
                num_mut = 0
                mutations = []
                while attempt < attempts and num_mut < self._mut_num and num_mut < iteration:
                    position = random.choice(tuple(PositionsSetUnion))
                    if position in mutations:
                        continue
                    old_gene = new_protein.genes[position - 1]

                    # no changes for CHARGED and TRP
                    set_name = 'Set1'
                    for name, sites in sets.items():
                        if old_gene.value in sites:
                            set_name = name
                            continue
                    pull = random.choices(probs[set_name][0], weights=probs[set_name][1][:-1])[0]
                    if random.random() > probs[set_name][1][-1]:
                        new_value = old_gene.value
                    else:
                        new_value = random.choice(pulls[pull])

                    new_gene = Gene(value=new_value)
                    new_protein.update_gene(position - 1, new_gene)

                    # проверка стабильности белка
                    if self.is_stable_protein(new_protein)\
                            and new_protein.num_changes <= iteration:
                        num_mut += 1
                        mutations.append(position)
                    else:
                        # Restore old gene
                        new_protein.update_gene(position - 1, old_gene)
                    attempt += 1
            new_population.append(new_protein)

        self._logger(f"Mutation: I will try to change {num_of_changed} proteins... {num_of_changed} proteins changed\n")

        self._population = new_population

    # crossover function
    def crossover(self, attempts=1):
        new_population = []
        for_cross = []  # белки для кроссовера

        for protein in self._population:
            new_protein = copy(protein)
            # условие на кроссинговер
            if random.random() < self._cros_prob:
                for_cross.append(new_protein)
            else:
                new_population.append(new_protein)

        # проверка на четное число белков в списке на кроссовер
        if len(for_cross) % 2 == 1:
            new_population.append(for_cross.pop())

        random.shuffle(for_cross)

        need = 0
        real = 0

        pair_cros_prob = 0.5  # crossover pair probability
        # цикл кроссовера(перемешивания генов белков)
        for protein1, protein2 in zip(for_cross[0:-1:2], for_cross[1::2]):
            need += 2
            new_protein1, new_protein2 = protein1, protein2
            for attempt in range(attempts):
                attempt_protein1, attempt_protein2 = copy(new_protein1), copy(new_protein2)
                mut_num = 0
                for i, (gene1, gene2) in enumerate(zip(attempt_protein1.genes, attempt_protein2.genes)):
                    if mut_num > self._mut_num:
                        continue
                    if random.random() < pair_cros_prob:
                        new_gene1 = Gene(value=gene2.value)
                        new_gene2 = Gene(value=gene1.value)
                        attempt_protein1.update_gene(i, new_gene1)
                        attempt_protein2.update_gene(i, new_gene2)
                        mut_num += 1

                if self.is_stable_protein(attempt_protein1) and self.is_stable_protein(attempt_protein2):
                    new_protein1 = attempt_protein1
                    new_protein2 = attempt_protein2
                    real += 2
                    break

            new_population.append(new_protein1)
            new_population.append(new_protein2)

        self._logger(f"Crossover: I will try to change {need} proteins... {real} proteins changed\n")

        self._population = new_population

    # selection function
    def selection(self, eval_param, save_n_best):
        def distribution(p, m):
            def evaluate(p: float, n: int) -> float:
                return p * pow(1 - p, n - 1)

            vs = []
            for j in range(1, m + 1):
                v = 0
                for j in range(1, j + 1):
                    v += evaluate(p, j)
                vs.append(v)
            return vs

        new_population = []
        pop_size = len(self._population)

        mean_value, mean_215_207 = self.mean_values(self._population)
        population = sorted(self._population, key=lambda x: (x.flatret, self.fitness(x, mean_value, mean_215_207)), reverse=True)  # increasing value

        for i in range(save_n_best):
            protein = copy(population[i])
            new_population.append(protein)

        q = distribution(eval_param, pop_size)
        for _ in range(pop_size - save_n_best):
            n, r = 0, random.uniform(0, q[-1])
            while r > q[n]:
                n += 1
            protein = copy(population[n])
            new_population.append(protein)

        new_population = sorted(new_population, key=lambda x: (x.flatret, self.fitness(x, mean_value, mean_215_207)))[0:pop_size]
        random.shuffle(new_population)

        self._population = new_population

    # функция подсчета value
    def compute(self):
        proteins_for_computing = []

        # Find existing calcs
        for protein in self._population:
            if protein.sequence not in self._computed:
                proteins_for_computing.append(protein)
        # Print to output file
        with open(".tempfile", "w") as ouf:
            for protein in proteins_for_computing:
                for idx, g1, g2 in protein.get_differences():
                    ouf.write(f"{g1}/{idx}/{g2} ")
                ouf.write("\n")
        os.rename(".tempfile", self._output_file)

        
        # Wait results
        while not os.path.exists(self._input_file):
            time.sleep(1)
        

        # Read results and save
        with open(self._input_file) as inf:
            value = []
            pka_215 = []
            pka_207 = []
            for line in inf.readlines():
                values = line.split()
                flat_ret = 1
                if len(values) <= 1:
                    value.append(50)
                    pka_215.append(-50)
                    pka_207.append(50)
                    flat_ret = 0
                else:
                    value.append(float(values[0]))
                    pka_215.append(float(values[1]))
                    pka_207.append(float(values[2]))
            for i, protein in enumerate(proteins_for_computing):
                self.save_computing(protein.sequence, value[i], pka_215[i], pka_207[i], flat_ret)

        # Write values to proteins
        for protein in self._population:
            values = self._computed[protein.sequence]
            protein.set_value(values[0])
            protein.set_pka215(values[1])
            protein.set_pka207(values[2])
            protein.set_flatret(values[3])

        # Remove out/inp filess
        os.remove(self._output_file)
        os.remove(self._input_file)

    def save_computing(self, sequence, value, pka_215, pka_207, flatret):
        if sequence not in self._computed:
            self._computed[sequence] = (value, pka_215, pka_207, flatret)
            with open(self._save_file, 'a') as f:
                f.write(f"{sequence} {value} {pka_215} {pka_207} {flatret}\n")

    # функция проверки стабильности белка(выполнения constraints)
    def is_stable_protein(self, protein):
        if self._checker is not None:
            return self._checker.check(protein)
        return True

    # функция выбора лучшего белка в популяции
    def get_best_protein(self):
        mean_value, mean_215_207 = self.mean_values(self._population)
        best_protein = max(self.population, key=lambda x: (x.flatret, self.fitness(x, mean_value, mean_215_207)))
        return best_protein, mean_value, mean_215_207

    # функция создания первой популяции
    def generate_population(self, default_sequence, default_value, default_pka215, default_pka207, pop_size,
                            from_computed=True):
        population = []

        self.save_computing(default_sequence, default_value, default_pka215, default_pka207, flatret=1)
        if from_computed:
            for sequence, value in self._computed.items():
                protein = Protein.create_protein(sequence, default_sequence, value=value[0],
                                                 pka215=value[1], pka207=value[2], flatret=value[3])
                population.append(protein)

            mean_value, mean_215_207 = self.mean_values(population)
            population = sorted(population, key=lambda x: (x.flatret, self.fitness(x, mean_value, mean_215_207)), reverse=True)[:pop_size]

        while len(population) < pop_size:
            protein = Protein.create_protein(default_sequence, default_sequence, value=default_value,
                                             pka215=default_pka215, pka207=default_pka207, flatret=1)
            population.append(protein)

        self._population = population

    def load_computed_proteins(self):
        if os.path.exists(self._save_file):
            with open(self._save_file, "r") as inf:
                for line in inf.readlines():
                    values = line.split()
                    self._computed[values[0]] = (float(values[1]), float(values[2]), float(values[3]), bool(values[4]))

    def print_current_population(self):
        for protein in self._population:
            self._logger(f"{protein.sequence}, BLA {protein.value} pka215 {protein.pka215}"
                         f" pka207 {protein.pka207} num of changes {protein.num_changes}\n")
            for idx, g1, g2 in protein.get_differences():
                self._logger(f"{g1}/{idx}/{g2} ")
            self._logger("\n")

    def fitness(self, prot, mean_value, mean_215_207):
        norm_value = 1
        delta = (prot.pka215-prot.pka207)/mean_215_207
        diff = abs(delta - norm_value)
        if diff != 0:
            fit = -prot.value/mean_value*self._weightval + 1/diff*self._weight215_207
        else:
            fit = -prot.value / mean_value * self._weightval + 1000 * self._weight215_207
        return fit

    def mean_values(self, population):
        mean_value = abs(sorted(population, key=lambda x: abs(x.value), reverse=True)[0].value)
        mean_prot = sorted(population, key=lambda x: abs(x.pka215 - x.pka207), reverse=True)[0]
        mean_215_207 = abs(mean_prot.pka215-mean_prot.pka207)
        return mean_value, mean_215_207

    def initial_step(self, population):
        ini_step = 1
        for protein in self._population:
            muts = len(protein.get_differences())
            if muts > ini_step:
                ini_step = muts
        return ini_step
