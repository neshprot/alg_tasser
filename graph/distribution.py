import matplotlib.pyplot as plt
import numpy as np
def distribution(p, m):
    def evaluate(p: float, n: int) -> float:
        return p * pow(1 - p, n - 1)

    vs = []
    for j in range(1, m + 1):
        v = 0
        for j in range(1, j + 1):
            v += evaluate(p, j)
        vs.append(1-v)
    return vs

q = distribution(0.5, 20)
x = np.linspace(1, 20, 20)

plt.plot(x, q)
plt.xlabel('номер индивида')
plt.ylabel('вероятность')
plt.xlim([1, 20])
plt.show()
