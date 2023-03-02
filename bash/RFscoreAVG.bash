#!/bin/bash
echo "This script is intended to average the docking energies obtained by RFScore over sevaral different poses and receptor conformations."
echo
echo "Please, run this script from the root directory of the ensemble-docking assay"
echo
echo "Please state the base directory name for the different receptor conformations."
echo "(The script assumes your folders names are composed of a base name and an integer number starting from 1)"
read basedirname
echo "Please state the base name of the ligands (e.g. ZINC)" 
read ligands

rootdir=`pwd`

for f in `ls -d ${basedirname}*/`; do
    cd ${f}
    ls -d ${ligands}*/ | sed -e s/"\/"// > ${rootdir}/ligandlist.out
    cd ${rootdir}
done 

########################################################################
# Average Average Binding Energies Over Different Receptor Conformations
########################################################################

for i in `cat ligandlist.out`; do
    echo ${i} >> names2.tmp
    counter=`cat ${basedirname}*/rfscorev3.out | grep ${i} | wc -l`
    cat ${basedirname}*/rfscorev3.out | grep ${i} | awk '{print $2}' | paste -sd+ | bc | awk -v var=$counter '{print $1/var}' >> average2.tmp
done

paste names2.tmp average2.tmp > EnsembleRFScores.tmp

##########################################
# Sort According to Average Binding Energy
##########################################

sort -r -k 2 EnsembleRFScores.tmp > EnsembleRFScores.out

rm names2.tmp average2.tmp EnsembleRFScores.tmp

exit


