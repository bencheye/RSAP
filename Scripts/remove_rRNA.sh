db_path=$1
core_num=$2
R1=$3
R2=$4
analysis_dir=$5
sample=$6
output1=$7
output2=$8
bowtie2 -x ${db_path} -p ${core_num}  --quiet -1 ${R1} -2 ${R2} --un-conc-gz ${analysis_dir}/${sample}
mv ${analysis_dir}/${sample}.1 ${output1}
mv ${analysis_dir}/${sample}.2 ${output2}
