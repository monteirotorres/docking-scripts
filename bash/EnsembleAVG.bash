#!/bin/bash
echo "This script is intended to average the docking energies obtained by Vina over sevaral different poses and receptor conformations."
echo
echo "Please, run this script from the root directory of the ensemble-docking assay"
echo
echo "Please state the base directory name for the different receptor conformations."
echo "(The script assumes your folders names are composed of a base name and an integer number starting from 1)"
read basedirname
echo "Please state the base name of the ligands (e.g. ZINC)" 
read ligands

rootdir=`pwd`

counter=0


###########################################
# Average Binding Energies of Several Poses
########################################### 

for f in `ls -d ${basedirname}*/`; do
    let counter=counter+1
    echo "Evaluating average binding energy for ligands in directory "${basedirname}${counter}"..."
    cd ${f}
    ls -d ${ligands}*/ | sed -e s/"\/"// > ${rootdir}/ligandlist.out

    for i in `ls -d ${ligands}*/`; do
        echo ${i} >> names.tmp
        numposes=`cat ${i}/log.txt | awk '/-----+/ {flag=1;next} /Writing/{flag=0} flag {print $2 }' | wc -l`
        cat ${i}/log.txt | awk '/-----+/ {flag=1;next} /Writing/{flag=0} flag {print $2 }' | paste -sd+ | bc | awk -v var=$numposes '{print $1/var}' >> average.tmp
        cat ${i}/log.txt | grep " 1 " | awk '{print $2}' >> best.tmp
        echo ${numposes} >> numposes.tmp
    done
    
    echo -e "Ligand\tAverage\tBest\tNºPoses" > scores.out
    paste names.tmp average.tmp best.tmp numposes.tmp >> scores.tmp
    sort -r -k 2 scores.tmp >> scores.out
    rm -rf names.tmp average.tmp best.tmp numposes.tmp scores.tmp
    cd ${rootdir}
done 

########################################################################
# Average Average Binding Energies Over Different Receptor Conformations
########################################################################

for i in `cat ligandlist.out`; do
    echo ${i} >> names2.tmp
    cat ${basedirname}*/scores.out | grep ${i} | awk '{print $2}' | paste -sd+ | bc | awk -v var=$counter '{print $1/var}' >> average2.tmp
done

paste names2.tmp average2.tmp > EnsembleScores.tmp

##########################################
# Sort According to Average Binding Energy
##########################################

sort -r -k 2 EnsembleScores.tmp > EnsembleScores.out

rm names2.tmp average2.tmp EnsembleScores.tmp

exit


