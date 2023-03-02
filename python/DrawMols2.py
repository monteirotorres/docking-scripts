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
                    help='Defines whether there is a partial input to be loaded\n')


args = parser.parse_args()

glutamine = "[CH2X4][CH2X4][CX3](=[OX1])[NX3H2]"
peptidomimetic = "[OX1]=CNC"
glutamine_mol = Chem.MolFromSmarts(glutamine)
peptidomimetic_mol = Chem.MolFromSmarts(peptidomimetic)
app = QApplication(sys.argv)
svgWidget = QSvgWidget()
for line in open(args.input_smiles):
    if not line.startswith(">"):
        smi = line.split()[1]
        id = line.split()[0]
        mol = Chem.MolFromSmiles(smi)
        # mh = Chem.AddHs(mol)
        # AllChem.EmbedMolecule(mh, AllChem.ETKDG())
        #
        # _, res = rdEHTTools.RunMol(mh)
        # static_chgs = res.GetAtomicCharges()[:mol.GetNumAtoms()]

        AllChem.Compute2DCoords(mol)

        gln_hits = list(mol.GetSubstructMatches(glutamine_mol))
        pep_hits = list(mol.GetSubstructMatches(peptidomimetic_mol))
        all_hits = set(chain.from_iterable(gln_hits + pep_hits))
        print(gln_hits)
        print(pep_hits)
        print(all_hits)

        gln_bonds = []
        for bond in glutamine_mol.GetBonds():
            for match in gln_hits:
                aid1 = match[bond.GetBeginAtomIdx()]
                aid2 = match[bond.GetEndAtomIdx()]
                gln_bonds.append(mol.GetBondBetweenAtoms(aid1,aid2).GetIdx())
        print(gln_bonds)

        pep_bonds = []
        for bond in peptidomimetic_mol.GetBonds():
            for match in pep_hits:
                aid1 = match[bond.GetBeginAtomIdx()]
                aid2 = match[bond.GetEndAtomIdx()]
                pep_bonds.append(mol.GetBondBetweenAtoms(aid1,aid2).GetIdx())
        print(pep_bonds)
        all_bonds = set(gln_bonds + pep_bonds)
        print(all_bonds)
        red = (0.8,0.5,0.5)
        blue = (0.5,0.5,0.8)
        lilac = (0.8,0.5,0.8)

        at_colors = {}
        for i in all_hits:
            if i in set(chain.from_iterable(gln_hits)) and i not in set(chain.from_iterable(pep_hits)):
                at_colors[i] = red
            if i in set(chain.from_iterable(pep_hits)) and i not in set(chain.from_iterable(gln_hits)):
                at_colors[i] = blue
            if i in set(chain.from_iterable(pep_hits)) and i in set(chain.from_iterable(gln_hits)):
                at_colors[i] = lilac
        print(at_colors)

        bd_colors = {}
        for i in all_bonds:
            if i in gln_bonds and i not in pep_bonds:
                bd_colors[i] = red
            if i in pep_bonds and i not in gln_bonds:
                bd_colors[i] = blue
            if i in pep_bonds and i in gln_bonds:
                bd_colors[i] = lilac
        print(bd_colors)


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

    if input("Show next mol?") == "n":
        break
    else:
        pass
