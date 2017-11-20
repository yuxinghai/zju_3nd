import os
work_dir="/data2/zhoulab/yuxinghai/zju/3nd_analysis/results"
anno_dir="/data2/zhoulab/yuxinghai/zju/anno/ce11"
bw_dir="/data2/zhoulab/yuxinghai/zju/3nd_analysis/results/bw"
raw_dir="/data2/zhoulab/yuxinghai/zju/3nd_data/ANCZJ170117_ZJ170117-02_2017-10-05_BHC57YCCXY/Cleandata"
script_dir="/data2/zhoulab/yuxinghai/zju/3nd_analysis/bin/script"
figure_dir="/data2/zhoulab/yuxinghai/zju/3nd_analysis/results/03_featurecount/figure"
results_dir="/data2/zhoulab/yuxinghai/zju/3nd_analysis/results/03_featurecount"
gene_SCF_dir="/data2/zhoulab/yuxinghai/zju/3nd_analysis/geneSCF-master-source-v1.1-p2"
samples = ["903","905","907","909","911","913","915","917","923","935","937","939","941","943","945","947","904","906","908","910","912","914","916","918","924","936","938","940","942","944","946","948"]
gene_list=["up","down","sig"]
exp_level=["up_Gene_sym","down_Gene_sym","sig_Gene_sym"]
period=["h_2","h_8","h_12","h_24"]
Tea_Dirctory_dir="/data2/zhoulab/yuxinghai/zju/anno/TEA_directory"
chrom_size="/data2/zhoulab/yuxinghai/zju/anno/ce11/ce11.chrom.sizes"


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
        expand("{figure_dir}/h_12_sample_distance.pdf",figure_dir=figure_dir),
        expand("{figure_dir}/h_8_sample_distance.pdf",figure_dir=figure_dir),
        expand("{figure_dir}/h_2_sample_distance.pdf",figure_dir=figure_dir),
	expand("{figure_dir}/h_24_sample_distance.pdf",figure_dir=figure_dir),
	expand("{results_dir}/gene_name/gene_SCF/h_2/down_Gene_sym",results_dir=results_dir),
	expand("{results_dir}/gene_name/gene_SCF/h_8/down_Gene_sym",results_dir=results_dir),
	expand("{results_dir}/gene_name/gene_SCF/h_12/down_Gene_sym",results_dir=results_dir),
	expand("{results_dir}/gene_name/gene_SCF/h_24/down_Gene_sym",results_dir=results_dir),
        expand("{results_dir}/DEgene/h_24/ano_downgene.csv",results_dir=results_dir),
        expand("{results_dir}/DEgene/h_24/ano_upgene.csv",results_dir=results_dir),
        expand("{results_dir}/DEgene/h_12/ano_downgene.csv",results_dir=results_dir),
        expand("{results_dir}/DEgene/h_12/ano_upgene.csv",results_dir=results_dir),
        expand("{results_dir}/DEgene/h_2/ano_downgene.csv",results_dir=results_dir),
        expand("{results_dir}/DEgene/h_2/ano_upgene.csv",results_dir=results_dir),
        expand("{results_dir}/DEgene/h_8/ano_downgene.csv",results_dir=results_dir),
        expand("{results_dir}/DEgene/h_8/ano_upgene.csv",results_dir=results_dir),
	expand("{work_dir}/02_mapping/sorted/{dataset}/{dataset}_uniq_sort.bam.bai",work_dir=work_dir,dataset=samples),
        expand("{work_dir}/bw/{dataset}.wig",work_dir=work_dir,dataset=samples),
        expand("{work_dir}/bw/{dataset}.bw",work_dir=work_dir,dataset=samples),
	expand("{results_dir}/gene_name/TEA/{period}/{gene_list}_tissue.svg",results_dir=results_dir,period=period,gene_list=gene_list),
	expand("{results_dir}/gene_name/TEA/{period}/{gene_list}_go.svg",results_dir=results_dir,period=period,gene_list=gene_list),
	expand("{results_dir}/gene_name/TEA/{period}/{gene_list}_phenotype.svg",results_dir=results_dir,period=period,gene_list=gene_list)







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
        "samtools sort -@ 20 {input} -o {output}"


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
        bam2wig.py -i {input} -o {bw_dir}/{wildcards.dataset} -t 1000000000 -s {chrom_size} -u
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
        "{figure_dir}/h_12_sample_distance.pdf",
        "{figure_dir}/h_2_sample_distance.pdf",
        "{figure_dir}/h_8_sample_distance.pdf",
        "{figure_dir}/h_24_sample_distance.pdf",
        
    shell:
        "Rscript {script_dir}/distance.R --latency-wait"


rule DEgene_filter_and_plot_h2:
    input:
        expand("{work_dir}/03_featurecount/{dataset}_featureCounts.txt",work_dir=work_dir, dataset=samples)

    output:
        "{results_dir}/figure/h_2/MAplot.png",
        "{results_dir}/figure/h_2/DE_number.pdf",
        "{results_dir}/gene_name/h_2/up",
        "{results_dir}/gene_name/h_2/down",
        "{results_dir}/gene_name/h_2/sig",
        "{results_dir}/gene_name/gene_SCF/h_2/up_Gene_sym",
        "{results_dir}/gene_name/gene_SCF/h_2/down_Gene_sym",
        "{results_dir}/DEgene/h_2/ano_downgene.csv",
        "{results_dir}/DEgene/h_2/ano_upgene.csv"

    shell:
        "Rscript {script_dir}/h2_deseq_refilter.R --latency-wait"

rule DEgene_filter_and_plot_h8:
    input:
        expand("{work_dir}/03_featurecount/{dataset}_featureCounts.txt",work_dir=work_dir, dataset=samples)

    output:
        "{results_dir}/figure/h_8/MAplot.png",
        "{results_dir}/figure/h_8/DE_number.pdf",
        "{results_dir}/gene_name/h_8/up",
        "{results_dir}/gene_name/h_8/down",
        "{results_dir}/gene_name/h_8/sig",
        "{results_dir}/gene_name/gene_SCF/h_8/up_Gene_sym",
        "{results_dir}/gene_name/gene_SCF/h_8/down_Gene_sym",
        "{results_dir}/DEgene/h_8/ano_downgene.csv",
        "{results_dir}/DEgene/h_8/ano_upgene.csv"

    shell:
        "Rscript {script_dir}/h8_deseq_refilter.R --latency-wait"

rule DEgene_filter_and_plot_h12:
    input:
        expand("{work_dir}/03_featurecount/{dataset}_featureCounts.txt",work_dir=work_dir, dataset=samples)

    output:
        "{results_dir}/figure/h_12/MAplot.png",
        "{results_dir}/figure/h_12/DE_number.pdf",
        "{results_dir}/gene_name/h_12/up",
        "{results_dir}/gene_name/h_12/down",
        "{results_dir}/gene_name/h_12/sig",
        "{results_dir}/gene_name/gene_SCF/h_12/up_Gene_sym",
        "{results_dir}/gene_name/gene_SCF/h_12/down_Gene_sym",
        "{results_dir}/DEgene/h_12/ano_downgene.csv",
        "{results_dir}/DEgene/h_12/ano_upgene.csv"

    shell:
        "Rscript {script_dir}/h12_deseq_refilter.R --latency-wait"

rule DEgene_filter_and_plot_h24:
    input:
        expand("{work_dir}/03_featurecount/{dataset}_featureCounts.txt",work_dir=work_dir, dataset=samples)

    output:
        "{results_dir}/figure/h_24/MAplot.png",
        "{results_dir}/figure/h_24/DE_number.pdf",
        "{results_dir}/gene_name/h_24/up",
        "{results_dir}/gene_name/h_24/down",
        "{results_dir}/gene_name/h_24/sig",
        "{results_dir}/gene_name/gene_SCF/h_24/up_Gene_sym",
        "{results_dir}/gene_name/gene_SCF/h_24/down_Gene_sym",
        "{results_dir}/DEgene/h_24/ano_downgene.csv",
        "{results_dir}/DEgene/h_24/ano_upgene.csv"

    shell:
        "Rscript {script_dir}/h24_deseq_refilter.R --latency-wait"





rule TEA:
    input:
        "{results_dir}/gene_name/{period}/{gene_list}"

    output:
        "{results_dir}/gene_name/TEA/{period}/{gene_list}_tissue.svg",
        "{results_dir}/gene_name/TEA/{period}/{gene_list}_go.svg",
        "{results_dir}/gene_name/TEA/{period}/{gene_list}_phenotype.svg"

    shell:
        """
        bash {script_dir}/tea.sh {wildcards.period}
        """
