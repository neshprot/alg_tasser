from collections import Counter

name = f'all_best'
words = ''
with open(name, 'r') as inf:
    s = 0
    counter = 0
    for i, line in enumerate(inf.readlines()):
        if len(line) > 1:
            first_word = line.split()[0]
            if first_word != 'Iteration:' and first_word != 'Crossover:' and first_word != 'Mutation:' and first_word != 'Step/Stop':
                if first_word == 'Current':
                    counter = i % 2
                    continue
                if counter == 0:
                    if i % 2 == 0:
                        s += 1
                        words += line
                else:
                    if i % 2 != 0:
                        s+=1
                        words += line

a =Counter(words.split())
sorted_a = sorted(a.items(), key=lambda x: x[1], reverse=True)
ox = []
oy = []

for value in sorted_a:
    ox.append(value[0])
    oy.append(float(value[1])/s*100)

with open('freq_and_repl', 'w') as file:
    for repl, freq in zip(ox, oy):
        file.write(f'{repl} {freq}\n')
