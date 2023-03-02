#!/bin/sh

echo "Script for organizing the vina results"
echo 
echo "You must have the vina_screen_get_top.py as an executable file on your path, or current directory"
echo
echo "How many molecules you want to colect?"
read molnum 

echo "Molecule" > names.tmp
echo "Best" > bestenergies.tmp 
echo "Worst" > worstenergies.tmp
rm -rf VinaReportTop*

for i in `vina_screen_get_top.py $molnum |  grep out.pdbqt`; do	 
	echo $i | sed s/"\/out.pdbqt"// >> names.tmp;
	head -2 $i | grep 'REMARK VINA RESULT' | awk '{print $4}' >> bestenergies.tmp;	
	sed '1!G;h;$!d' $i | awk '1;/MODEL/{exit}' | grep 'REMARK VINA RESULT' | awk '{print $4}' >> worstenergies.tmp
	awk '1;/ENDMDL/{exit}' $i >> VinaReportTop$molnum.pdbqt;
done

paste names.tmp bestenergies.tmp worstenergies.tmp > VinaReportTop$molnum.txt
rm -rf names.tmp *energies.tmp

exit
