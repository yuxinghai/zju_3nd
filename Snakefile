import os
work_dir="/data2/zhoulab/yuxinghai/zju/3nd_analysis/results"
anno_dir="/data2/zhoulab/yuxinghai/zju/anno/ce11"
bw_dir="/data2/zhoulab/yuxinghai/zju/3nd_analysis/results/bw"
raw_dir="/data2/zhoulab/yuxinghai/zju/3nd_data/ANCZJ170117_ZJ170117-02_2017-10-05_BHC57YCCXY/Cleandata"
script_dir="/data2/zhoulab/yuxinghai/zju/3nd_analysis/bin/script"
figure_dir="/data2/zhoulab/yuxinghai/zju/3nd_analysis/results/03_featurecount/figure"
results_dir="/data2/zhoulab/yuxinghai/zju/3nd_analysis/results/03_featurecount"
gene_SCF_dir="/data2/zhoulab/yuxinghai/zju/3nd_analysis/geneSCF-master-source-v1.1-p2"
samples = ["911","912","913","914","915","916","917","918"]
gene_list=["up","down","sig"]
exp_level=["up_Gene_sym","down_Gene_sym","sig_Gene_sym"]
Tea_Dirctory_dir="/data2/zhoulab/yuxinghai/zju/anno/TEA_directory"


rule all:
    input:
        expand("{work_dir}/01_fastqc/{dataset}/{dataset}_R1_fastqc.html",work_dir=work_dir, dataset=samples),
        expand("{work_dir}/01_fastqc/{dataset}/{dataset}_R2_fastqc.html",work_dir=work_dir, dataset=samples),
        expand("{anno_dir}/rRNA/SAindex",anno_dir=anno_dir),
        expand("{anno_dir}/genome/SAindex",anno_dir=anno_dir),
        expand("{work_dir}/02_mapping/rRNA/{dataset}/{dataset}Unmapped.out.mate1", work_dir=work_dir, dataset=samples),
        expand("{work_dir}/02_mapping/rRNA/{dataset}/{dataset}Unmapped.out.mate2", work_dir=work_dir, dataset=samples),
        expand("{work_dir}/02_mapping/rRNA/{dataset}/{dataset}Unmapped.out.mate1_sort", work_dir=work_dir, dataset=samples),
        expand("{work_dir}/02_mapping/rRNA/{dataset}/{dataset}Unmapped.out.mate2_sort", work_dir=work_dir, dataset=samples),
	expand("{work_dir}/02_mapping/sorted/{dataset}/{dataset}_SortAligned.out.sam",work_dir=work_dir, dataset=samples),
	expand("{work_dir}/02_mapping/sorted/{dataset}/{dataset}_uniq.bam",work_dir=work_dir, dataset=samples),
	expand("{work_dir}/02_mapping/sorted/{dataset}/{dataset}_uniq_sort.bam",work_dir=work_dir, dataset=samples),
        expand("{work_dir}/03_featurecount/{dataset}_featureCounts.txt", work_dir=work_dir, dataset=samples),
        expand("{work_dir}/03_featurecount/{dataset}_featureStat.log", work_dir=work_dir, dataset=samples),
        expand("{figure_dir}/sample_distance.pdf", figure_dir=figure_dir),
        expand("{figure_dir}/damage_vs_control_MAplot.png", figure_dir=figure_dir),
        expand("{figure_dir}/Gene_number_in_DE_level.pdf", figure_dir=figure_dir),
        expand("{results_dir}/gene_name/up", results_dir=results_dir),
        expand("{results_dir}/gene_name/down", results_dir=results_dir),
        expand("{results_dir}/gene_name/sig", results_dir=results_dir),
        expand("{results_dir}/gene_name/gene_SCF/up_Gene_sym", results_dir=results_dir),
        expand("{results_dir}/gene_name/gene_SCF/down_Gene_sym", results_dir=results_dir),
        expand("{results_dir}/DEgene/damage_vs_control_downgene.csv", results_dir=results_dir),
        expand("{results_dir}/DEgene/damage_vs_control_upgene.csv", results_dir=results_dir),
        expand("{results_dir}/DEgene/ano_damage_vs_control_downgene.csv", results_dir=results_dir),
        expand("{results_dir}/DEgene/ano_damage_vs_control_upgene.csv", results_dir=results_dir),
        expand("{results_dir}/gene_name/gene_SCF/sig_Gene_sym", results_dir=results_dir),
        expand("{results_dir}/gene_name/gene_SCF/{exp_level}", results_dir=results_dir,exp_level=exp_level),
        expand("{results_dir}/gene_name/gene_SCF_result/{exp_level}/{exp_level}_GO_all_wb_functional_classification.tsv", results_dir=results_dir,exp_level=exp_level),
        expand("{results_dir}/gene_name/gene_SCF_result/{exp_level}/{exp_level}_KEGG_cel_enrichment_plot.pdf", results_dir=results_dir,exp_level=exp_level),
        expand("{results_dir}/gene_name/TEA/{gene_list}_tissue.svg", results_dir=results_dir,gene_list=gene_list),
        expand("{results_dir}/gene_name/TEA/{gene_list}_go.svg", results_dir=results_dir,gene_list=gene_list),
        expand("{results_dir}/gene_name/TEA/{gene_list}_phenotype.svg", results_dir=results_dir,gene_list=gene_list)


rule fastqc:
    input:
        R1 = raw_dir+"/{dataset}/{dataset}_R1.fq.gz",
        R2 = raw_dir+"/{dataset}/{dataset}_R2.fq.gz",

    output:
        "{work_dir}/01_fastqc/{dataset}/{dataset}_R1_fastqc.html",
        "{work_dir}/01_fastqc/{dataset}/{dataset}_R2_fastqc.html",
        "{work_dir}/01_fastqc/{dataset}/{dataset}_R1_fastqc.zip",
        "{work_dir}/01_fastqc/{dataset}/{dataset}_R2_fastqc.zip"

    shell:
        "fastqc -f fastq -t 2 --noextract -o {work_dir}/01_fastqc/{wildcards.dataset} {input.R1} {input.R2}"

rule build_index:
    input:
        rRNA_fa="{anno_dir}/ce11_rRNA.fa",
        gtf="{anno_dir}/c_elegans.PRJNA13758.WS258.canonical_geneset.gtf",
        fa="{anno_dir}/ce11.fa"

    output:
        rRNA_index="{anno_dir}/rRNA/SAindex",
        genome_index="{anno_dir}/genome/SAindex"

    threads: 20

    shell: """
        STAR --runThreadN 20 --runMode genomeGenerate --genomeDir {anno_dir}/rRNA --genomeFastaFiles {input.rRNA_fa} --genomeSAindexNbases 6 --outFileNamePrefix {anno_dir}/rRNA/ce11_rRNA

        STAR --runThreadN 20 --runMode genomeGenerate --genomeDir {anno_dir}/genome --genomeFastaFiles {input.fa} --sjdbGTFfile {input.gtf} --sjdbOverhang 149 --genomeSAindexNbases 12 --outFileNamePrefix {anno_dir}/genome/ce11_ws258
        """

rule rRNA_map:
    input:
        R1 = raw_dir+"/{dataset}/{dataset}_R1.fq.gz",
        R2 = raw_dir+"/{dataset}/{dataset}_R2.fq.gz"

    output:
        "{work_dir}/02_mapping/rRNA/{dataset}/{dataset}Unmapped.out.mate1",
        "{work_dir}/02_mapping/rRNA/{dataset}/{dataset}Unmapped.out.mate2"

    threads: 20

    shell:
        """STAR --genomeDir {anno_dir}/rRNA --readFilesIn {input.R1} {input.R2} --readFilesCommand gunzip -c --runThreadN 20 --outFileNamePrefix {work_dir}/02_mapping/rRNA/{wildcards.dataset}/{wildcards.dataset} --outReadsUnmapped Fastx --outFilterMatchNmin 40 --outFilterScoreMinOverLread 0 --outFilterMatchNminOverLread 0
	"""


rule sort_STAR_result:
    input:
        un1="{work_dir}/02_mapping/rRNA/{dataset}/{dataset}Unmapped.out.mate1",
        un2="{work_dir}/02_mapping/rRNA/{dataset}/{dataset}Unmapped.out.mate2"


    output:
        fastq1="{work_dir}/02_mapping/rRNA/{dataset}/{dataset}Unmapped.out.mate1_sort",
        fastq2="{work_dir}/02_mapping/rRNA/{dataset}/{dataset}Unmapped.out.mate2_sort"

    shell:
        """
        cat {input.un1}| paste -d " " - - - - | sort -k1,1 -S 10G | tr ' ' '\n' > {output.fastq1}
        cat {input.un2}| paste -d " " - - - - | sort -k1,1 -S 10G | tr ' ' '\n' > {output.fastq2}
        """

rule genome_map:
    input:
        f_fastq1="{work_dir}/02_mapping/rRNA/{dataset}/{dataset}Unmapped.out.mate1_sort",
        f_fastq2="{work_dir}/02_mapping/rRNA/{dataset}/{dataset}Unmapped.out.mate2_sort"


    output:
        f_Sam="{work_dir}/02_mapping/sorted/{dataset}/{dataset}_SortAligned.out.sam"

    threads: 20

    shell:
        """
        STAR --genomeDir {anno_dir}/genome --readFilesIn {input.f_fastq1} {input.f_fastq2} --runThreadN 20 --outFileNamePrefix {work_dir}/02_mapping/sorted/{wildcards.dataset}/{wildcards.dataset}_Sort --outReadsUnmapped Fastx --outFilterMatchNmin 40 --outFilterScoreMinOverLread 0 --outFilterMatchNminOverLread 0
        """

rule sam_2_uniquebam:
    input:
        "{work_dir}/02_mapping/sorted/{dataset}/{dataset}_SortAligned.out.sam"

    output:
        u_bam="{work_dir}/02_mapping/sorted/{dataset}/{dataset}_uniq.bam"

    shell:
        "samtools view -@ 8 -h -bq 255 {input} > {output.u_bam}"

rule uniquebam_sort_and_index:
    input:
        "{work_dir}/02_mapping/sorted/{dataset}/{dataset}_uniq.bam"

    output:
        "{work_dir}/02_mapping/sorted/{dataset}/{dataset}_uniq_sort.bam"

    threads: 20

    shell:
        "samtools view -@ 20 {input} -o {output}"


rule bam_index_and_2bw:
    input:
        "{work_dir}/02_mapping/sorted/{dataset}/{dataset}_uniq_sort.bam"

    output:
        "{work_dir}/02_mapping/sorted/{dataset}/{dataset}_uniq_sort.bam.bai",
        "{work_dir}/bw/{dataset}.wig",
        "{work_dir}/bw/{dataset}.bw"

    shell:
        """
        samtools index {input}
        bam2wig.py -i {input} -o {bw_dir}/{dataset} -t 1000000000 -s {chrom_size} -u
        """

rule feature_count:
    input:
        "{work_dir}/02_mapping/sorted/{dataset}/{dataset}_uniq_sort.bam"

    output:
        "{work_dir}/03_featurecount/{dataset}_featureCounts.txt",
        "{work_dir}/03_featurecount/{dataset}_featureStat.log"

    shell:
        "Rscript {script_dir}/featurecount.R"

rule sample_distance_plot:
    input:
        expand("{work_dir}/03_featurecount/{dataset}_featureCounts.txt",work_dir=work_dir, dataset=samples)

    output:
        "{figure_dir}/sample_distance.pdf"

    shell:
        "Rscript {script_dir}/distance.R --latency-wait"


rule DEgene_filter_and_plot:
    input:
        expand("{work_dir}/03_featurecount/{dataset}_featureCounts.txt",work_dir=work_dir, dataset=samples)

    output:
        "{results_dir}/figure/damage_vs_control_MAplot.png",
        "{results_dir}/figure/Gene_number_in_DE_level.pdf",
        "{results_dir}/gene_name/up",
        "{results_dir}/gene_name/down",
        "{results_dir}/gene_name/sig",
        "{results_dir}/gene_name/gene_SCF/up_Gene_sym",
        "{results_dir}/gene_name/gene_SCF/down_Gene_sym",
        "{results_dir}/DEgene/damage_vs_control_downgene.csv",
        "{results_dir}/DEgene/damage_vs_control_upgene.csv",
        "{results_dir}/DEgene/ano_damage_vs_control_downgene.csv",
        "{results_dir}/DEgene/ano_damage_vs_control_upgene.csv"

    shell:
        "Rscript {script_dir}/damage_deseq_refilter.R"

rule prepare_geneSCF_input:
    input:
        up="{results_dir}/gene_name/gene_SCF/up_Gene_sym",
        down="{results_dir}/gene_name/gene_SCF/down_Gene_sym"

    output:
        sig="{results_dir}/gene_name/gene_SCF/sig_Gene_sym"
    shell:
        """
        cat {input.up} {input.down} >{output.sig}
        """


rule gene_SCF:
    input:
        "{results_dir}/gene_name/gene_SCF/{exp_level}"

    output:
        "{results_dir}/gene_name/gene_SCF_result/{exp_level}/{exp_level}_GO_all_wb_functional_classification.tsv",
        "{results_dir}/gene_name/gene_SCF_result/{exp_level}/{exp_level}_KEGG_cel_enrichment_plot.pdf"


    shell:
        "{gene_SCF_dir}/geneSCF_batch"


rule TEA:
    input:
        "{results_dir}/gene_name/{gene_list}"

    output:
        "{results_dir}/gene_name/TEA/{gene_list}_tissue.svg",
        "{results_dir}/gene_name/TEA/{gene_list}_go.svg",
        "{results_dir}/gene_name/TEA/{gene_list}_phenotype.svg"

    shell:
        """
        bash {script_dir}/tea.sh
        mv down_* {results_dir}/gene_name/TEA
        mv up_* {results_dir}/gene_name/TEA
        mv sig_* {results_dir}/gene_name/TEA
        """
