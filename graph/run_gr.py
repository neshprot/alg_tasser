import matplotlib.pyplot as plt

run_num = '6'
inp_file = f'data3_{run_num}'
probs = [None, None]
config_weight = [None, None]
ox = []
BLA = []
pka215 = []
pka207 = []

with open(inp_file) as file:
    lines = file.readlines()
    params = lines[0].split()
    probs[0] = params[0]
    probs[1] = params[1]
    config_weight[0] = params[2]
    config_weight[1] = params[3]
    for i, line in enumerate(lines[1::]):
        values = line.split()
        ox.append(i+1)
        BLA.append(float(values[0]))
        pka215.append(float(values[1]))
        pka207.append(float(values[2]))

delta_pka = [pka215[i] - pka207[i] for i in range(len(pka215))]

fig, (ax1, ax2) = plt.subplots(1, 2)
ax1.set_title(f'cros prob: {probs[0]}, mut prob: {probs[1]}')
ax2.set_xlabel(f'population, weight for pka: {config_weight[1]}')
ax2.set_ylabel('pka values')
ax2.plot(ox, pka215, label='pka215')
ax2.plot(ox, pka207, label='pka207')
ax2.plot(ox, delta_pka, 'ko', label='pka215 - pka207')

ax1.set_xlabel(f'population, weight for BLA: {config_weight[0]}')
ax1.set_ylabel('BLA')
ax1.bar(ox, BLA, label='BLA')
fig.set_figheight(5)
fig.set_figwidth(10)
plt.legend()
plt.savefig(f'data3_graph{run_num}.png')
