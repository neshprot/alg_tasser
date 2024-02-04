import matplotlib.pyplot as plt

mutations = {}

with open('start_set', 'r') as inf:
    for line in inf.readlines():
        key = line.split('/')[1]
        values = line.split()
        mutations[key] = mutations.get(key, []) + [(values[3], values[0])]

print(mutations)
ox = []
oy = []
text = []
for num in mutations.keys():
    for value in mutations[num]:
        ox.append(int(num))
        oy.append(float(value[0]))
        text.append(value[1])

for i in range(len(ox)):
    if ox[i] != 60 and ox[i] != 196:
        plt.annotate(text[i], (ox[i], oy[i]+0.04))
    else:
        print(f'plt.annotate(\'{text[i]}\', ({ox[i]}, {oy[i]}+0.04))')

plt.annotate('T/196/G', (196, -7.619999999999999+0.01))
plt.annotate('T/196/A', (196, -7.68))
plt.annotate('T/196/M', (196, -7.6899999999999995-0.05))
plt.annotate('T/196/Q', (196, -7.6899999999999995-0.12))
plt.annotate('T/196/C', (196, -7.6899999999999995-0.18))
plt.annotate('T/196/L', (196, -7.699999999999999-0.23))
plt.annotate('T/196/I', (196, -7.720000000000001-0.3))
plt.annotate('T/196/V', (196, -7.73-0.37))

plt.annotate('Y/60/F', (60, -8.799999999999999+0.04))
plt.annotate('Y/60/M', (60, -8.97-0.06))
plt.annotate('Y/60/L', (60, -8.989999999999998-0.1))
plt.annotate('Y/60/Q', (60, -9.01+0.04))

plt.scatter(ox, oy, s=5)
plt.show()
