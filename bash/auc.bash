#!/bin/bash
for i in `seq 1 ${points}`;do
    ax=`sed -n ${i}p a.out | awk '{print $1}'`
    ay=`sed -n ${i}p a.out | awk '{print $2}'`
    let i=i+1
    bx=`sed -n ${i}p a.out | awk '{print $1}'`
    by=`sed -n ${i}p a.out | awk '{print $2}'`
    echo ${ax} ${ay} ${bx} ${by}
done > aucpoints.tmp

auc=`awk '{printf "%11.10f\n", (($3-$1)*($2+$4))/2}' aucpoints.tmp | paste -sd+ | bc`

echo "The area under the curve is "${auc}
