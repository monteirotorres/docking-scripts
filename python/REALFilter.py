import os
import multiprocess as mp
global filters
global real_db

min_mw = 240
max_mw = 480
min_slogp = 0
max_slogp = 3
min_hba = 2
max_hba = 10
min_hbd = 2
max_hbd = 5
min_rotb = 3
max_rotb = 8
min_fsp3 = 0.1
max_fsp3 = 0.5
min_tpsa = None
max_tpsa = None
min_qed = None
max_qed = None
pains = True
brenk = True
zinc = True
lilly = True
nih = True
type = 'S'

filters = {'mw': (min_mw, max_mw), 'slogp': (min_slogp, max_slogp),
           'hba': (min_hba, max_hba), 'hbd': (min_hbd, max_hbd),
           'rotb': (min_rotb, max_rotb), 'fsp3': (min_fsp3, max_fsp3),
           'psa': (min_tpsa, max_tpsa), 'qed': (min_qed, max_qed),
           'pains': pains, 'brenk': brenk, 'zinc': zinc,
           'lilly': lilly, 'nih': nih, 'type': type}

rm_filters = []
for filter, value in filters.items():
    if value is None or value == (None, None):
        rm_filters.append(filter)

for filter in rm_filters:
    del filters[filter]

real_db = "/home/torres/work/Molecules/REAL"


def filter_line(line):
    vars = {}
    vars['smiles'], vars['idnumber'], vars['reagent1'], vars['reagent2'], vars['reagent3'], vars['reagent4'], vars['reaction'],	vars['mw'],	vars['slogp'], vars['hba'], vars['hbd'], vars['rotbonds'], vars['fsp3'], vars['tpsa'], vars['qed'], vars['pains'], vars['brenk'], vars['nih'], vars['zinc'], vars['lilly'], vars['lead_like'], vars['lead_like_350'], vars['fragments'], vars['strict_fragments'], vars['ppi_modulators'], vars['natural_product_like'], vars['type'] = line.split('\t')
    for var, value in vars.items():
        if var in filters:
            if var == 'type':
                if value == filters[var]:
                    continue
                else:
                    return
            elif var in ['pains', 'brenk', 'lilly', 'zinc', 'nih']:
                if bool(value):
                    return
                else:
                    continue
            else:
                if float(filters[var][0]) < float(value) < float(filters[var][1]):
                    continue
                else:
                    return
    return line

def filter_file(file):
    with open(os.path.join(real_db, file),'r') as source:
        next(source)
        for line in source:
            if filter_line(line.strip()):
                print(line, flush=True)

filelist = []
for file in os.listdir(real_db):
    if file.endswith('smiles'):
        filelist.append(file)


p = mp.Pool(12)
for filtered_file in p.imap_unordered(filter_file, filelist):
    pass
