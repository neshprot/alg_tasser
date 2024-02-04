import os
import time
from logger import FileLogger

logger = FileLogger('logout')

def compute(sets):
    # Wait results
    with open(".tempfile", "w") as ouf:
        for mutation in sets:
            ouf.write(f'{mutation} ')
            ouf.write("\n")
    os.rename(".tempfile", 'tmp_out')

    while not os.path.exists('tmp_inp'):
        time.sleep(1)

    with open('tmp_inp', 'r') as inf:
        i = 0
        for line in inf.readlines():
            if len(line) == 0:
                logger(f'{sets[i]} None\n')
            else:
                values = line.split()
                pka_215 = float(values[0])
                pka_207 = float(values[1])
                logger(f'{sets[i]} pka215 {pka_215} pka207 {pka_207} pka215-pka207 {pka_215-pka_207}\n')
            i += 1
    os.remove('tmp_inp')
    os.remove('tmp_out')


all_muts = []

with open('main_tmp') as ouf:
    for line in ouf.readlines():
        all_muts.append(line.split()[0])

sets = []
for i, mut in enumerate(all_muts):
    sets.append(mut)
    if (i+1) % 20 == 0:
        compute(sets)
        sets = []
