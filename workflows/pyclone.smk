import pandas as pd

samples = pd.read_csv("samples/samples_info_curated.txt",sep="\t").set_index("sample", drop=False)

rule all:
    input:
        inter_samples = expand("data/pyclone_output/{sample}", sample =samples.index)


rule intersect:
    input:
        snv = "data/snv/{sample}_SNV.tsv.gz",
        cnv = "data/cnv/{sample}_CNV.tsv.gz"
    output:
        out_dir = "data/intersect/{sample}_intersect.txt"

    shell:
        "bedtools intersect -a {input.snv} -b {input.cnv} -wa -wb > {output.out_dir}"

rule pyclone_input:
    input:
        "data/intersect/{sample}_intersect.txt"
    output:
        "data/pyclone_input/{sample}_pyclone.txt"
    shell:
        "Rscript scripts/get_pyclone_input_script.R {input} {output}"

rule run_pyclone:
    input:
       "data/pyclone_input/{sample}_pyclone.txt"
    output:
       "data/pyclone_output/{sample}"
    run:
        os.mkdir(output[0])
	pur = samples.loc[wildcards.sample,'purity']
	shell("PyClone run_analysis_pipeline --in_files {input} --working_dir {output} --tumour_contents {pur}  --density pyclone_binomial --min_cluster_size 2")
