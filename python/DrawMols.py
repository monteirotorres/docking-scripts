#!/usr/bin/env python
from rdkit.Chem import AllChem
from rdkit import Chem
from rdkit.Chem.Draw import rdMolDraw2D
import sys
from PyQt5.QtWidgets import QApplication
import argparse
from PyQt5.QtSvg import QSvgWidget
from rdkit.Chem import rdEHTTools
from rdkit.Chem.Draw import SimilarityMaps
import time
from itertools import chain

parser = argparse.ArgumentParser()

parser.add_argument('-f',
                    dest='input_smiles',
                    action='store',
                    default=None,
                    help='Smiles file from which to draw mols\n')
parser.add_argument('-id',
                    dest='mol_id',
                    type=str,
                    default=None,
                    help='Target mol to draw\n')


args = parser.parse_args()

#query1 = "[CH2X4][CH2X4][CX3](=[OX1])[NX3H2]"
#query2 = "[OX1]=CNC"
# query1 = "[#7;A;H2X3][#6;X3](=[O;X1])[#6;A;H2X4][#6;A;X4][#6]-[#6]=O"
# #query1 = "[NX3;H2,H1][#6;X3](=[O;X1])[#6][#6;!R;A;H2X4][#6]-[#6]=O"
# query2 = '[#6!R]-[#6!R]-[#6R]-1-[#6R]-[#6R]-[#7HR]-[#6R]-1=O'


query1 = "[#7;A;H2X3][#6;X3](=[O;X1])[*!R][*!R][*!R]-[*]=O" # GLN Analog
query2 = '[*!R]-[*!R]-[*R]-1-[*R]-[*R]-[#7HR]-[#6R]-1=O'# 5 Ring GLN Analog
query3 = '[*!R]-[*!R]-[*R]-1-[*R]-[*R]-[*R]-[#7HR]-[#6R]-1=O'# 6 Ring GLN Analog

query1_mol = Chem.MolFromSmarts(query1)
if query2:
    query2_mol = Chem.MolFromSmarts(query2)
if query3:
    query3_mol = Chem.MolFromSmarts(query3)
app = QApplication(sys.argv)
svgWidget = QSvgWidget()
for line in open(args.input_smiles):
    if not line.startswith(">"):
        smi = line.split()[1]
        id = line.split()[0]
        if args.mol_id:
            if args.mol_id == id:
                pass
            else:
                continue
        print('\n'+id)
        print(smi+'\n')
        mol = Chem.MolFromSmiles(smi)
        # mh = Chem.AddHs(mol)
        # AllChem.EmbedMolecule(mh, AllChem.ETKDG())
        #
        # _, res = rdEHTTools.RunMol(mh)
        # static_chgs = res.GetAtomicCharges()[:mol.GetNumAtoms()]

        AllChem.Compute2DCoords(mol)

        query1_hits = list(mol.GetSubstructMatches(query1_mol))
        if query2:
            query2_hits = list(mol.GetSubstructMatches(query2_mol))
        else:
            query2_hits = []
        if query3:
            query3_hits = list(mol.GetSubstructMatches(query3_mol))
        else:
            query3_hits = []
        all_hits = set(chain.from_iterable(query1_hits + query2_hits + query3_hits))
        # print(query1_hits)
        # print(query2_hits)
        # print(all_hits)

        query1_bonds = []
        for bond in query1_mol.GetBonds():
            for match in query1_hits:
                aid1 = match[bond.GetBeginAtomIdx()]
                aid2 = match[bond.GetEndAtomIdx()]
                query1_bonds.append(mol.GetBondBetweenAtoms(aid1,aid2).GetIdx())
        # print(query1_bonds)

        query2_bonds = []
        if query2:
            for bond in query2_mol.GetBonds():
                for match in query2_hits:
                    aid1 = match[bond.GetBeginAtomIdx()]
                    aid2 = match[bond.GetEndAtomIdx()]
                    query2_bonds.append(mol.GetBondBetweenAtoms(aid1,aid2).GetIdx())
        # print(query2_bonds)

        query3_bonds = []
        if query3:
            for bond in query3_mol.GetBonds():
                for match in query3_hits:
                    aid1 = match[bond.GetBeginAtomIdx()]
                    aid2 = match[bond.GetEndAtomIdx()]
                    query3_bonds.append(mol.GetBondBetweenAtoms(aid1,aid2).GetIdx())
        all_bonds = set(query1_bonds + query2_bonds + query3_bonds)
        # print(all_bonds)
        red = (0.8,0.5,0.5)
        blue = (0.5,0.5,0.8)
        green = (0.5,0.8,0.5)
        rb = (0.8,0.5,0.8)
        rg = (0.8,0.8,0.5)
        gb = (0.5,0.8,0.8)
        rgb = (0.8,0.8,0.8)


        at_colors = {}
        for i in all_hits:
            if i in set(chain.from_iterable(query1_hits)) and i not in set(chain.from_iterable(query2_hits+query3_hits)):
                at_colors[i] = red
            elif i in set(chain.from_iterable(query2_hits)) and i not in set(chain.from_iterable(query1_hits+query3_hits)):
                at_colors[i] = green
            elif i in set(chain.from_iterable(query3_hits)) and i not in set(chain.from_iterable(query1_hits+query2_hits)):
                at_colors[i] = blue
            elif i in set(chain.from_iterable(query1_hits)) and i in set(chain.from_iterable(query2_hits)) and i not in set(chain.from_iterable(query3_hits)):
                at_colors[i] = rg
            elif i in set(chain.from_iterable(query1_hits)) and i in set(chain.from_iterable(query3_hits)) and i not in set(chain.from_iterable(query2_hits)):
                at_colors[i] = rb
            elif i in set(chain.from_iterable(query2_hits)) and i in set(chain.from_iterable(query3_hits)) and i not in set(chain.from_iterable(query1_hits)):
                at_colors[i] = gb
            elif i in set(chain.from_iterable(query1_hits)) and i in set(chain.from_iterable(query2_hits)) and i in set(chain.from_iterable(query3_hits)):
                at_colors[i] = rgb
        # print(at_colors)

        bd_colors = {}
        for i in all_bonds:
            if i in query1_bonds and i not in query2_bonds and i not in query3_bonds:
                bd_colors[i] = red
            elif i in query2_bonds and i not in query1_bonds and i not in query3_bonds:
                bd_colors[i] = green
            elif i in query3_bonds and i not in query1_bonds and i not in query2_bonds:
                bd_colors[i] = blue
            elif i in query1_bonds and i in query2_bonds and i not in query3_bonds:
                bd_colors[i] = rg
            elif i in query1_bonds and i in query3_bonds and i not in query2_bonds:
                bd_colors[i] = rb
            elif i in query2_bonds and i in query3_bonds and i not in query1_bonds:
                bd_colors[i] = gb
            elif i in query1_bonds and i in query2_bonds and i in query3_bonds:
                bd_colors[i] = rgb

        # print(bd_colors)


        drawer = rdMolDraw2D.MolDraw2DSVG(800, 600)
        drawer.DrawMolecule(mol, highlightAtoms=all_hits,
                            highlightAtomColors=at_colors,
                            highlightBonds=all_bonds,
                            highlightBondColors=bd_colors)
        #SimilarityMaps.GetSimilarityMapFromWeights(mol, list(static_chgs), draw2d=drawer)
        drawer.FinishDrawing()
        svg = drawer.GetDrawingText().replace('svg:', '')


        svg_bytes = bytearray(svg, encoding='utf-8')

        svgWidget.renderer().load(svg_bytes)
        svgWidget.show()
    else:
        continue
    if args.mol_id:
        input("Press ENTER to exit...")
        break
    if input("Show next mol?") == "n":
        break
    else:
        pass
