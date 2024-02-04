from collections import Counter

import matplotlib.pyplot as plt

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

a = set([int(i.split('/')[1]) for i in words.split()])
with open('replacements', 'w') as file:
    for site in sorted(a):
        file.write(f'{str(site)} ')
