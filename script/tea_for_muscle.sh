#!/usr/bin/bash
des_dir=/data2/zhoulab/yuxinghai/zju/3nd_analysis/results/03_featurecount/gene_name/$1
dic_dir=/data2/zhoulab/yuxinghai/zju/anno/TEA_directory
out_dir=/data2/zhoulab/yuxinghai/zju/3nd_analysis/results/03_featurecount/gene_name/TEA/$1
NAMES=(transition_down transition_up  ms_siggene ms_not_siggene) 
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
