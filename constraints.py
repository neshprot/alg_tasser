from functools import partial

import numpy as np
from numpy.linalg import norm

from data import Protein


class Constraints:
    """
    Класс, соедржащий констрэинты
    """

    def __init__(self):
        self._constraints = []
        self._minv = 0
        self._maxv = 10

    def add(self, f) -> None:
        """
        Добавляет новый констрэинт.
        Все функции должна обладать сигнатурой T -> bool.

        :param f: констрэинт.
        :return:
        """
        self._constraints.append(f)

    def check(self, x) -> bool:
        """
        Проверяет удовлетворение объекта всем констрэинтам

        :param x:
        :return: True, если все функции ограничения вернули True
        """
        for f in self._constraints:
            if not f(x):
                return False
        return True

    @property
    def minv(self):
        return self._minv
    
    @property
    def maxv(self):
        return self._maxv,

    @minv.setter
    def minx(self, x):
        self._minv = x
        
    @maxv.setter
    def maxv(self, x):
        self._maxv = x


def constraint_max_charge(p: Protein, max_charge: float) -> bool:
    return abs(p.charge) < max_charge


def constraint_n_charged(p: Protein, max_n_charged: float) -> bool:
    n = len(p.genes)
    count = 0
    for i in range(0, n):
        if p[i].charged:
            count += 1
    return count < max_n_charged


def constraint_distances(p: Protein, min_distance: float, coords: np.ndarray, positions_set) -> bool:
    assert len(p.genes) == len(coords)
    n = len(p.genes)
    for i in range(0, n):
        for j in range(i + 1, n):
            if (i + 1) in positions_set and (j + 1) in positions_set:
                if p.genes[i].charged and p.genes[j].charged:
                    dist = norm(coords[i] - coords[j])
                    if dist < min_distance:
                        return False
    return True


def constraint_pka(p: Protein, position_set, value: float):
    minav = p.minv
    maxav = p.maxv
    for idx, g1, g2 in p.get_differences():
        if idx == 215:
            value = p.value
            if (isinstance(value, float)):
                if value <= 5 and value > minav:
                    minav = value
                if value >= 5 and value < maxav:
                    maxav = value
                if value < p.minv or value > p.maxv:
                    p.minv = minav
                    p.maxv = maxav
                    return False
    p.minv = minav
    p.maxv = maxav
    return True


def constraint_included(p: Protein, aminoacids_set, positions_set) -> bool:
    for n, gene in enumerate(p.genes):
        if n + 1 in positions_set:
            if gene.value in aminoacids_set:
                return True
    return False


def constraint_max_num_changes(p: Protein, max_num_changes) -> bool:
    return p.num_changes <= max_num_changes


if __name__ == "__main__":
    from evolution import PositionsSetUnion, PositionsSet2
    import random
    from data import AMINOACIDS

    constraint = Constraints()

    f1 = partial(constraint_distances, min_distance=5.0, coords=coordinates, positions_set=PositionsSetUnion)
    f2 = partial(constraint_max_charge, max_charge=7)
    f3 = partial(constraint_n_charged, max_n_charged=60)
    f4 = partial(constraint_pka, position_set=PositionsSet2, value = 0)

    constraint.add(f1)
    constraint.add(f2)
    constraint.add(f3)
    constraint.add(f4)

    import random
    from data import AMINOACIDS

    sequence = random.choices(AMINOACIDS, k=200)
    protein = Protein(sequence)

    results_f = [f(protein) for f in [f1, f2, f3]]
    print(f"Constraint 1: {results_f[0]}")
    print(f"Constraint 2: {results_f[1]}")
    print(f"Constraint 3: {results_f[2]}")
    print(f"Constraint 4: {results_f[3]}")

    results_checker = constraint.check(protein)
    print(f"Constraints checker: {results_checker}")

    assert all(results_f) == results_checker

    print("OK!")
