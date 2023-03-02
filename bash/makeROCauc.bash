#!/bin/bash

echo "##############################################"
echo "#              Make ROC + AUC                #"
echo "#   Creates a table containing Sensitivity   #"
echo "#  and 1-Specificity for each tested cutoff  #"
echo "#       TORRES, P.H.M. and GOMES, D.E.       #"
echo "#                26-11-2017                  #"
echo "##############################################"

export LC_NUMERIC=en_US.UTF-8

echo "
Please provide a file in which the first column correspond to the name of the ligand and the second column corresponds to the binding energy"
read -e -p "mainfile: " mainfile

echo "
Please state the fix part of the actives name"
read -e -p "Actives string: " actives

echo "
Please state the fix part of the decoys name"
read -e -p "Decoys string: " decoys

echo "
How many points should there be in your ROC curve?"
read -e -p "Number of points: " points

echo "
Do you want to generate a xvg file?"
read -e -p "Generate xvg?: " xvg

let points=${points}-1

# Sort by binding energy
sort -n -k2 ${mainfile} > mainsorted.tmp

#### Write active file
fgrep ${actives} ${mainfile} > actives.tmp
fgrep ${decoys}  ${mainfile} > decoys.tmp


# Best energy, worst energy, range and step
best=`head -n1 mainsorted.tmp | awk '{print $2}'`
worst=`tail -n1 mainsorted.tmp | awk '{print $2}'`
range=`echo|awk -v best=$best -v worst=$worst '{print best-worst}' OFMT='%.6f'`
step=`echo|awk -v points=$points -v range=$range '{print range/points}' OFMT='%.6f'`

#Number of active and decoys
activetot=`wc -l  actives.tmp |awk '{print $1}'`
decoytot=` wc -l  decoys.tmp  |awk '{print $1}'`

echo -e "
     best = ${best}
    worst =  ${worst}
    range = ${range}
     step = ${step}
activetot = ${activetot}
 decoytot = ${decoytot}
"

cutoff=$worst

let points=${points}+1

for i in `seq 1 $points`; do
	tn=0
	tp=0
	echo ${cutoff} >> cutoff.tmp
	
     a=0
    a=`awk -v cutoff=${cutoff} '$2<=cutoff' actives.tmp |wc -l`
    a=`echo ${a}`
    let tp=${tp}+${a}

     d=0
    d=`awk -v cutoff=${cutoff} '$2>cutoff' decoys.tmp |wc -l`
    d=`echo ${d}`
    let tn=${tn}+${d}

  echo ${cutoff} ${tp} ${tn} 
        
  echo "$tp/$activetot" | bc -l >> sensitivity.tmp
    echo "1-($tn/$decoytot)" | bc -l >> 1-specificity.tmp
	cutoff=`echo $worst + $i*$step | bc`
done
echo -e "Cutoff\t1-Specificity\tSensitivity" > ROC.txt
paste cutoff.tmp 1-specificity.tmp sensitivity.tmp >> ROC.txt

tac ROC.txt > ROCinv.tmp

if [ ${xvg} == "yes" -o ${xvg}=="y" ] ; then
echo -e "
@    title \"ROC Curve\"
@    xaxis  label \"1-Specificity\"
@    yaxis  label \"Sensitivity\"
@TYPE xy
@    s0 symbol 2
@    s0 symbol size 0.3
@    s0 linestyle 1
" > ROC.xvg
paste 1-specificity.tmp sensitivity.tmp >> ROC.xvg
fi

for i in `seq 1 ${points}`;do
    ax=`sed -n ${i}p ROCinv.tmp | awk '{print $1}'`
    ay=`sed -n ${i}p ROCinv.tmp | awk '{print $2}'`
    let i=i+1
    bx=`sed -n ${i}p ROCinv.tmp | awk '{print $1}'`
    by=`sed -n ${i}p ROCinv.tmp | awk '{print $2}'`
    echo ${ax} ${ay} ${bx} ${by}
done > aucpoints.tmp

auc=`awk '{printf "%11.10f\n", (($3-$1)*($2+$4))/2}' aucpoints.tmp | paste -sd+ | bc`

echo "For "${points}" points, The area under the curve is "${auc}"."

echo "Npoints"\t"AUC"\n${points}\t${auc} > auc.log

rm -rf sensitivity.tmp 1-specificity.tmp cutoff.tmp mainsorted.tmp ROCinv.tmp aucpoints.tmp

