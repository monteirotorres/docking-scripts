#!/usr/bin/env python

import os

workdir = '/home/torres/work/vitor/world'
file = os.path.join(workdir, 'world.mol2')

current_mol_block = []
first_mol = True
getmolname = False
known_mols = {}
with open(file, 'r') as fin:
    for line in fin.read().splitlines():
        if line.startswith('@<TRIPOS>MOLECULE'):
            if first_mol is True:
                first_mol = False
                pass
            else:
                current_mol = current_mol_block[1]
                if current_mol not in known_mols:
                    known_mols[current_mol] = 1
                else:
                    known_mols[current_mol] += 1
                outfile_name = current_mol+'_'+str(known_mols[current_mol])+'.mol2'
                print('\n\nWriting '+outfile_name+'\n')
                #print('\n'.join(current_mol_block))
                with open(os.path.join(workdir, outfile_name), 'w') as fout:
                    fout.writelines('\n'.join(current_mol_block))
            current_mol_block = []
            current_mol_block.append(line)
        else:
            current_mol_block.append(line)
