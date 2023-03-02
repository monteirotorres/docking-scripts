#!/bin/bash
echo "This script is intended to normalize the docking energy of a virtual screening run by dividing these energies by the cube root of the number of non-hydrogen atoms, as suggested by Pan, et al. (J Chem. Inf. Comput. Sci.2003,43,267-272)"
echo
echo "In order for this script to porpperly work you must have the Vega software installed and in your PATH."
echo 
echo "Please indicate a directory where your ligands are located"
read liganddir
echo "Please provide a file containing two columns in which the first one corresponds to the ligand name (pdbqt format) and the second corresponds to the docking energy. For Example:"
echo "ZINC06716957	-8.6"
echo "ZINC01612996	-8.59571"
echo "ZINC03978005	-8.50071"
echo "ZINC52955754	-8.42429"
echo "."
echo "."
echo "."
read energyfile
echo "Do you want the final results to be ordered according to the normalized energies?"
read order

while read line; do 
    ligandname=`echo "$line" | awk '{print $1}'` ;
    heavyatoms=`vega $liganddir/$ligandname.pdbqt -f INFO | grep -a Heavy | awk '{print $3}'` ;
    echo "$line" | awk -v var=$heavyatoms '{print $2/var^(1/3)}' >> normalized.tmp ;
done < $energyfile
paste <(awk '{print $1}' $energyfile) <(awk '{printf "%2.3f\n", $2}' $energyfile) <(awk '{printf "%2.3f\n", $1}' normalized.tmp) > NormalizedScores.tmp

if [ $order = yes -o $order = y ]; then
    cat NormalizedScores.tmp | sort -r -k 3 >NormalizedScores.out
else
    mv NormalizedScores.tmp NormalizedScores.out
fi


rm -rf normalized.tmp NormalizedScores.tmp
