#/bin/bash
inputData=$1
outputData=$2

echo -e "Geneid\tgene_type\tgene_name" > ${outputData}
less ${inputData} | grep -v '^##' | awk -F 'gene_id ' '{print $2}' | awk -F ';' '{print $1,$3,$4}' | grep -v ' level ' | uniq | sed -e 's/gene_type//g' -e 's/gene_name//g' -e 's/"//g' >> ${outputData}
