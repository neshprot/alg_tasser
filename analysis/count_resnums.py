import matplotlib.pyplot as plt
resids = {}
with open('first_resids', 'r') as file:
    for line in file.readlines():
        if int(line) in resids.keys():
            resids[int(line)] += 1
        else:
            resids[int(line)] = 1

ox = []
oy = []

for key, value in resids.items():
    ox.append(key)
    oy.append(value)
    print(f'{key} - {value}')
fig, ax = plt.subplots()
ax.bar(ox, oy)
ax.set_xlabel('res ids')
ax.set_ylabel('frequency')
plt.show()
