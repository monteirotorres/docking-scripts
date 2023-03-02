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
import inspect

real_db = "/home/torres/work/Molecules/REAL/*.smi"
smile_files = [filepath for filepath in glob.glob(real_db, recursive = True)]
query="[#7;A;HX3][#6;X3](=[O;X1])[#6;A;X4][#6;A;H2X4][#6]-[#6]=O" # Generic Glutamine Analog
#query2 = "O=[#6]-[#6]-[#6]-[#6]-1-[#6]-[#6]-[#7]-[#6]-1=O"

query_mol = Chem.MolFromSmarts(query)
#query_mol2 = Chem.MolFromSmarts(query2)


def analyse_smiles(line):
    smi = line.split()[0]
    id = line.split()[1]
    try:
        trial_mol = Chem.MolFromSmiles(smi)
        if trial_mol.HasSubstructMatch(query_mol): #or trial_mol.HasSubstructMatch(query_mol2):
            print(smi+'\t'+id, flush=True)
            return True
    except:
        pass
    return False


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
            for response in p.imap_unordered(analyse_smiles, source, chunksize=100):
                pass
end_time = datetime.now()
runtime = end_time - start_time
print('TOTAL RUNTIME: '+str(runtime)+' ms')
