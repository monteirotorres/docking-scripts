#! /bin/bash
echo "Script for running vina"
echo 
echo "Please state the fix part of the pdbqt file names"
read fixed
 
for f in $fixed*.pdbqt; do
    b=`basename $f .pdbqt`
    echo Processing ligand $b
    mkdir -p $b
    vina --config conf.txt --ligand $f --out ${b}/out.pdbqt --log ${b}/log.txt
done
