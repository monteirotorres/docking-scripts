#!/usr/bin/env python

import os
import multiprocess
import argparse
import glob
from rdkit import Chem
from rdkit.Chem.Fraggle import FraggleSim
from rdkit import DataStructs
from rdkit import RDLogger
from datetime import datetime
global glutamine
global query_mol

anodynes_db = "/home/torres/work/ZincSmiles/**/*.smi"
smile_files = [filepath for filepath in glob.glob(anodynes_db, recursive = True)]
grp = "[NH2][CH]=O"
grp_mol = Chem.MolFromSmiles(grp)
test = "CCC[C@H]([NH2+][C@@H](C)C(=O)Nc1cc(C)c(Cl)cc1OC)C([O-])=O"
test_mol = Chem.MolFromSmiles(test)
glutamine = "N[C@@H](CCC(N)=O)C(O)=O"
query_mol = Chem.MolFromSmiles(glutamine)
query_fp = Chem.RDKFingerprint(query_mol)
query_mol.HasSubstructMatch(grp_mol, )
test_mol.HasSubstructMatch(grp_mol)

def fragment_similarity(smiles):
    trial_mol = Chem.MolFromSmiles(smiles)
    fragments = FraggleSim.generate_fraggle_fragmentation(trial_mol)
    set_fragments = set()
    for fragment in fragments:
        for frag in fragment.split('.'):
            set_fragments.add(frag)

    try:
        for fragment in set_fragments:
            trial_fragment = Chem.MolFromSmarts(fragment)
            response = DataStructs.cDataStructs.SokalSimilarity(Chem.RDKFingerprint(trial_fragment), query_fp)
            if response > .3:
                return response
    except:
        pass
    return False


def analyse_smiles(line):
    smi = line.split()[0]
    id = line.split()[1]
    response = fragment_similarity(smi)
    if response:
        print(id+'\t'+smi+'\t'+str(round(response, 3)), flush=True)


parser = argparse.ArgumentParser()

parser.add_argument('-pi', '--partial-input',
                    dest='partial_input',
                    action='store',
                    default=None,
                    help='Defines whether there is a partial input to be loaded\n')

args = parser.parse_args()
skip_files = []
if args.partial_input:
    with open(args.partial_input, 'r') as fin:
        for line in fin:
            if line.startswith('>'):
                processed_file = line.split()[1]
                skip_files.append(processed_file)
    with open("previous.out", "w") as previous_out:
        with open(args.partial_input, 'r') as fin:
            for line in fin:
                if line.startswith('>') and skip_files[-1] in line:
                    break
                previous_out.write(line)

RDLogger.DisableLog('rdApp.*')
start_time = datetime.now()
p = multiprocess.Pool()
for smile_file in smile_files:
    if smile_file not in skip_files[:-1]:
        print("> "+smile_file, flush=True)
        with open(smile_file) as source:
            next(source)
            for response in p.imap_unordered(analyse_smiles, source):
                pass
end_time = datetime.now()
runtime = end_time - start_time
print('TOTAL RUNTIME: '+str(runtime)+' ms')
