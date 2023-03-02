#!/usr/bin/env python
# Imports
###############################################################################
import argparse
import textwrap as tw
from rdkit import Chem
from rdkit.Chem import AllChem
from multiprocessing import cpu_count
from rdkit.Chem.SaltRemover import SaltRemover, InputFormat
import multiprocess

# Description
###############################################################################
print(tw.dedent("""\

        #############################################################

                             ConvertSmiles
                  Copyright (C) 2018  Pedro H. M. Torres
                       The University of Cambridge

             This program comes with ABSOLUTELY NO WARRANTY

        Given a csv file containing a series of smiles, followed by
        their names, convert them into 3D entities, remove salt ions,
        optimize their structures and output in separated PDB or mol2
        files. Also outputs a sdf file if requested.

        ***requires rdkit

        #############################################################
        \n
      """))

# Global Variables
###############################################################################
parser = argparse.ArgumentParser()

parser.add_argument( dest='smiles_list',
                    type=str,
                    metavar='smiles.tsv',
                    help='csv file containing a series of smiles, followed by their names.\n')

parser.add_argument('-f', '--format',
                    dest='format',
                    type=str, default=None,
                    metavar='',
                    help='Output format [pdb or mol2].\n')

parser.add_argument('-s', '--sdf',
                    dest='sdf',
                    action='count',
                    default=0,
                    help='Generate sdf file.\n')

parser.add_argument('-np',
                    dest='np',
                    type=int,
                    default=cpu_count(),
                    help='Generate sdf file.\n')

global args

args = parser.parse_args()

remover = SaltRemover()


def analyse_line(line):
    try:
        smile = line.split('\t')
        m = smile[1]
        name = smile[0]
        print(name)
        m2 = Chem.MolFromSmiles(m)
        m3 = remover.StripMol(m2)
        m4 = Chem.AddHs(m3)
        AllChem.EmbedMolecule(m4,AllChem.ETKDG())
        AllChem.UFFOptimizeMolecule(m4)
        m4.SetProp("_Name",name)
        if args.format == 'pdb':
            w = Chem.PDBWriter(name+'.pdb')
            w.write(m4)
        elif args.format == 'mol2':
            w = Chem.SDWriter(name+'.mol2')
            w.write(m4)
        return m4
    except KeyboardInterrupt:
        exit()
    except:
        print('Failed to convert molecule'+str(smile))
        return None


# Main Function
###############################################################################
def main():
    p = multiprocess.Pool(args.np)
    mols = []
    source = open(args.smiles_list, 'r')
    for m in p.map(analyse_line, source):
        try:
            if args.sdf == 1:
                if m is not None:
                    mols.append(m)
        except:
            prnint("HERE")
    if args.sdf == 1:
        w = Chem.SDWriter('all.sdf')
        for m in mols: w.write(m)



# Execute
###############################################################################
main()
