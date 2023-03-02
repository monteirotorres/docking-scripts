#!/usr/bin/env python
import os
import multiprocess
import argparse
import glob
from rdkit import Chem
from rdkit import RDLogger
from datetime import datetime
global query_mol

smile_file = '/home/torres/work/Molecules/REAL-Filtered/Filtered_noblanks.smi'


# query="[NX3;H2,H1][#6;X3](=[O;X1])[#6][#6;!R;A;H2X4][#6]-[#6]=O" # Generic Glutamine Analog
query1 = "[#7;A;H2X3][#6;X3](=[O;X1])[*!R][*!R][*!R]-[*]=O" # GLN Analog
query2 = '[*!R]-[*!R]-[*R]-1-[*R]-[*R]-[#7HR]-[#6R]-1=O'# 5 Ring GLN Analog
query3 = '[*!R]-[*!R]-[*R]-1-[*R]-[*R]-[*R]-[#7HR]-[#6R]-1=O'# 6 Ring GLN Analog
#query = 'O=[#6]-[#6]-[#6]-[#6]-1-[#6]-[#6]-[#7]-[#6]-1=O'
#query = "[OX1]=CNC" # Peptidomimetic

query1_mol = Chem.MolFromSmarts(query1)
query2_mol = Chem.MolFromSmarts(query2)
query3_mol = Chem.MolFromSmarts(query3)

def analyse_smiles(line):
    smi = line.split()[0]
    id = line.split()[1]
    try:
        trial_mol = Chem.MolFromSmiles(smi)
        if trial_mol.HasSubstructMatch(query1_mol) or trial_mol.HasSubstructMatch(query2_mol) or trial_mol.HasSubstructMatch(query3_mol):
            print(smi+'\t'+id, flush=True)
            return True
    except:
        pass
    return False


RDLogger.DisableLog('rdApp.*')
start_time = datetime.now()
p = multiprocess.Pool()

with open(smile_file) as source:
    for response in p.imap_unordered(analyse_smiles, source, chunksize=100):
        pass
end_time = datetime.now()
runtime = end_time - start_time
print('TOTAL RUNTIME: '+str(runtime)+' ms')
