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

anodynes_db = "/home/torres/work/Molecules/Zinc200-500_0-3_270M/**/*.smi"
smile_files = [filepath for filepath in glob.glob(anodynes_db, recursive = True)]
#query="[#7;A;HX3][#6;X3](=[O;X1])[#6;A;X4][#6;A;H2X4][#6]-[#6]=O" # Generic Glutamine Analog
# query = '[#6!R]-[#6!R]-[#6R]-1-[#6R]-[#6R]-[#7HR]-[#6R]-1=O'
# query2 = '[H]N[C@@H](CC1CCN([H])C1=O)C=O'
query1 = "[#7;A;H2X3][#6;X3](=[O;X1])[*!R][*!R][*!R]-[*]=O" # GLN Analog
query2 = '[*!R]-[*!R]-[*R]-1-[*R]-[*R]-[#7HR]-[#6R]-1=O'# 5 Ring GLN Analog
query3 = '[*!R]-[*!R]-[*R]-1-[*R]-[*R]-[*R]-[#7HR]-[#6R]-1=O'# 6 Ring GLN Analog
query1_mol = Chem.MolFromSmarts(query1)
query2_mol = Chem.MolFromSmarts(query2)
query3_mol = Chem.MolFromSmarts(query3)



def analyse_smiles(line):
    smi = line.split()[0]
    id = line.split()[1]
    try:
        trial_mol = Chem.MolFromSmiles(smi)
        if trial_mol.HasSubstructMatch(query1_mol) or trial_mol.HasSubstructMatch(query3_mol) or trial_mol.HasSubstructMatch(query3_mol):
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
