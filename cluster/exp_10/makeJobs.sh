#!/bin/bash
while IFS= read -r line; do
    seed=`echo $line | awk '{print $1}'`
    L=`echo $line | awk '{print $2}'`
    ite=`echo $line | awk '{print $3}'`
    echo  $seed $L $ite
    sed -e "s|@seed@|$seed|g" oar_script.sh > ${L}_${ite}_submit.sh
    sed -i "s|@L@|$L|g" ${L}_${ite}_submit.sh
    sed -i "s|@ite@|$ite|g" ${L}_${ite}_submit.sh
    chmod +x ./${L}_${ite}_submit.sh
    oarsub -S ./${L}_${ite}_submit.sh
done < params.dat

