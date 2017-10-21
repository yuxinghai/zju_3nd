#!/usr/bin/bash
des_dir=/data2/zhoulab/yuxinghai/zju/3nd_analysis/results/03_featurecount/gene_name/
dic_dir=/data2/zhoulab/yuxinghai/zju/anno/TEA_directory
out_dir=/data2/zhoulab/yuxinghai/zju/DEgene/gene_name/TEA
NAMES=(down sig up) 
numLanes=${#NAMES[@]}
cd ${out_dir}
for (( i=0; i<${numLanes}; i++ )); do
    name=${NAMES[$i]}
    gene_list=${des_dir}/${name}
    tea -d ${dic_dir}/anatomy_dict_95_33_WS258.csv -s ${gene_list} ${name}_tissue tissue
done

for (( i=0; i<${numLanes}; i++ )); do
     name=${NAMES[$i]}
     gene_list=${des_dir}/${name}
     tea -d ${dic_dir}/go_dict_95_100_WS258.csv -s ${gene_list} ${name}_go go
done

for (( i=0; i<${numLanes}; i++ )); do
     name=${NAMES[$i]}
     gene_list=${des_dir}/${name}
     tea -d ${dic_dir}/phenotype_dict_95_50_WS258.csv -s ${gene_list} ${name}_phenotype phenotype
done
