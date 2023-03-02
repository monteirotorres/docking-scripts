#!/usr/bin/env python
# Imports
###############################################################################

from rdkit import Chem

# Main Function
###############################################################################
def main():
    inf = '/home/torres/work/coronavirus/docking/NS5-GLN-Analogs-REAL/glide/sp.sdf'
    suppl = Chem.SDMolSupplier(inf)
    for mol in suppl:
        if mol is None: continue

        w = Chem.PDBWriter(mol.GetProp('s_lp_Variant')+'.pdb')
        w.write(mol)



# Execute
###############################################################################
main()
