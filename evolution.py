import os
import random
import time
from abc import abstractmethod, ABC
from copy import copy

from data import Protein, Gene

PullAPlus = "STNQCWYEDH"
PullBPlus = "PGAVILMFEDH"
PullA = "STNQCWY"
PullB = "PGAVILMF"
PullC = "STNQCGPAVILMFYW"
PullD = "RKHED"

ResiduesSet1 = "STNQCGPAVILMFYWEDH"
ResiduesSet2 = "STNQCGPAVILMFYW"
ResiduesSet3 = "RHKDESTNQCGPAVILMFYW"

PositionsSet1 = {23, 56, 88, 89, 92, 93, 96, 188, 215}
PositionsSet2 = {52, 60, 86, 121, 122, 124, 125, 126, 128, 129, 140, 141, 142, 144, 145, 146, 147, 148, 181, 184, 185,
                 186, 189, 191, 192, 193, 211, 218}
PositionsSet3 = {18, 19, 20, 22, 24, 26, 27, 30, 48, 49, 50, 51, 53, 54, 55, 57, 58, 59, 61, 62, 63, 80, 81, 82, 84, 85,
                 87, 90, 91, 94, 95, 97, 98, 99, 100, 114, 115, 116, 117, 118, 119, 120, 127, 130, 131, 134, 136, 137,
                 138, 139, 143, 149, 150, 151, 155, 177, 178, 179, 180, 182, 183, 187, 190, 194, 195, 197, 198, 199,
                 207, 208, 209, 210, 212, 213, 214, 216, 217, 220, 221, 222, 224, 225, 226}
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


class ProteinEvolution(Evolution, BaseFunction):
    def __init__(self, population, mut_prob, mut_num, cros_prob, input_file, output_file, save_file, logger, checker=None):
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

    def mutation(self, attempts=1):
        """

        :param attempts: число попыток инциниализации на один protein
        :return:
        """

        new_population = []
        num_of_changed = 0
        first_p = 0.6
        second_p = 0.4

        for protein in self._population:
            new_protein = copy(protein)
            if random.random() < self._mut_prob:
                num_of_changed += 1
                attempt = 0
                num_mut = 0
                while num_mut < self._mut_num and attempt < attempts:
                    position = random.choice(tuple(PositionsSetUnion))
                    old_gene = new_protein.genes[position - 1]
                    new_value = old_gene.value

                    if position in PositionsSet1:
                        if old_gene.polared:
                            if random.random() < first_p:
                                new_value = random.choice(PullBPlus)
                            elif random.random() < second_p:
                                new_value = random.choice(PullAPlus)
                        else:
                            if random.random() < first_p:
                                new_value = random.choice(PullAPlus)
                            elif random.random() < second_p:
                                new_value = random.choice(PullBPlus)
                    elif position in PositionsSet2:
                        if old_gene.polared:
                            if random.random() < first_p:
                                new_value = random.choice(PullBPlus)
                            elif random.random() < second_p:
                                new_value = random.choice(PullA)
                        else:
                            if random.random() < first_p:
                                new_value = random.choice(PullAPlus)
                            elif random.random() < second_p:
                                new_value = random.choice(PullB)
                    elif position in PositionsSet3:
                        if old_gene.charged:
                            if random.random() < first_p:
                                new_value = random.choice(PullC)
                            elif random.random() < second_p:
                                new_value = random.choice(PullD)
                        else:
                            if random.random() < first_p:
                                new_value = random.choice(PullD)
                            elif random.random() < second_p:
                                new_value = random.choice(PullC)

                    new_gene = Gene(value=new_value)
                    new_protein.update_gene(position - 1, new_gene)

                    if self.is_stable_protein(new_protein):
                        num_mut += 1
                    else:
                        # Restore old gene
                        new_protein.update_gene(position - 1, old_gene)
                    attempt += 1
            new_population.append(new_protein)

        self._logger(f"Mutation: I will try to change {num_of_changed} proteins... {num_of_changed} proteins changed\n")

        self._population = new_population

    def crossover(self, attempts=1):
        new_population = []
        for_cross = []

        for protein in self._population:
            new_protein = copy(protein)
            if random.random() < self._cros_prob:
                for_cross.append(new_protein)
            else:
                new_population.append(new_protein)

        if len(for_cross) % 2 == 1:
            new_population.append(for_cross.pop())

        random.shuffle(for_cross)

        need = 0
        real = 0

        pair_cros_prob = 0.5  # crossover pair probability
        for protein1, protein2 in zip(for_cross[0:-1:2], for_cross[1::2]):
            need += 2
            new_protein1, new_protein2 = protein1, protein2
            for attempt in range(attempts):
                attempt_protein1, attempt_protein2 = copy(new_protein1), copy(new_protein2)
                for i, (gene1, gene2) in enumerate(zip(attempt_protein1.genes, attempt_protein2.genes)):
                    if random.random() < pair_cros_prob:
                        new_gene1 = Gene(value=gene2.value)
                        new_gene2 = Gene(value=gene1.value)
                        attempt_protein1.update_gene(i, new_gene1)
                        attempt_protein2.update_gene(i, new_gene2)

                if self.is_stable_protein(attempt_protein1) and self.is_stable_protein(attempt_protein2):
                    new_protein1 = attempt_protein1
                    new_protein2 = attempt_protein2
                    real += 2
                    break

            new_population.append(new_protein1)
            new_population.append(new_protein2)

        self._logger(f"Crossover: I will try to change {need} proteins... {real} proteins changed\n")

        self._population = new_population

    # TODO: расширить возможности этапа селекции
    def selection(self, eval_param, save_n_best):
        def distribution(p, m):
            def evaluate(p: float, n: int) -> float:
                return p * pow(1 - p, n - 1)

            vs = []
            for i in range(1, m + 1):
                v = 0
                for j in range(1, i + 1):
                    v += evaluate(p, j)
                vs.append(v)
            return vs

        new_population = []
        pop_size = len(self._population)

        population = sorted(self._population, key=lambda x: x.value, reverse=True)

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

        new_population = sorted(new_population, key=lambda x: x.value)[0:pop_size]
        random.shuffle(new_population)

        self._population = new_population

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

            # Fake computing
            # from run_simulate_computing import run_simulate_computing
            # run_simulate_computing()

        # Read results and save
        with open(self._input_file) as inf:
            for protein in proteins_for_computing:
                value = float(inf.readline())
                self.save_computing(protein.sequence, value)

        # Write values to proteins
        for protein in self._population:
            value = self._computed[protein.sequence]
            protein.set_value(value)

        # Remove out/inp filess
        os.remove(self._output_file)
        os.remove(self._input_file)

    def save_computing(self, sequence, value):
        if sequence not in self._computed:
            self._computed[sequence] = value
            with open(self._save_file, 'a') as f:
                f.write(f"{sequence} {value}\n")

    def is_stable_protein(self, protein):
        if self._checker is not None:
            return self._checker.check(protein)
        return True

    def get_best_protein(self) -> Protein:
        best_protein = max(self.population, key=lambda x: x.value)
        return best_protein

    def generate_population(self, default_sequence, default_value, pop_size, from_computed=True):
        population = []

        self.save_computing(default_sequence, default_value)

        if from_computed:
            for sequence, value in self._computed.items():
                protein = Protein.create_protein(sequence, default_sequence, value=value)
                population.append(protein)
            population = sorted(population, key=lambda x: x.value, reverse=True)[:pop_size]

        while len(population) < pop_size:
            protein = Protein.create_protein(default_sequence, default_sequence, value=default_value)
            population.append(protein)

        self._population = population

    def load_computed_proteins(self):
        if os.path.exists(self._save_file):
            with open(self._save_file, "r") as inf:
                for line in inf.readlines():
                    sequence, value = line.split()
                    self._computed[sequence] = float(value)

    def print_current_population(self):
        for protein in self._population:
            self._logger(f"{protein.sequence}, {protein.value}, {protein.num_changes}\n")
            for idx, g1, g2 in protein.get_differences():
                self._logger(f"{g1}/{idx}/{g2} ")
            self._logger("\n")
