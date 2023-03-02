#!/bin/bash
echo "Please provide a file in which the first column corresponds to the molecule ZINC code."
read zincnames
echo "How many of the top molecules should have their names retrieved from ZINC database?"
read nmolecules

for i in `head -${nmolecules} $zincnames | awk '{print$1}'` ; do 
    curl -s http://zinc.docking.org/substance/$i | grep " <abbr title=" | awk -F "\"" '{print $2}';
done
exit
