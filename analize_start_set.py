muts = []

with open('start_set', 'r') as inf:
    for line in inf.readlines():
        values = line.split()
        if len(values) == 2:
            continue
        muts.append((values[0], values[2], values[4], float(values[6])))
muts = list(set(muts))
new_set = sorted(muts, key=lambda x: x[3], reverse=True)
with open('final_set', 'w') as ouf:
    for values in new_set[:20:]:
        ouf.write(f'{values[0]} {values[1]} {values[2]} {values[3]}\n')
