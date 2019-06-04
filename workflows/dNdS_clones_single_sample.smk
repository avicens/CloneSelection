import pandas as pd

samples = pd.read_csv("samples_prueba_info.txt",sep="\t").set_index("sample", drop=False)

rule all:
    input:
        #inter_samples = expand("data/pyclone_output/{sample}", sample =samples.index),
        #input_loci = expand("data/pyclone_input/{sample}_pyclone.txt", sample = samples.index),
        general_output = expand("data/raxml_trees/{sample}_raxml.log", sample = samples.index)

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
        "data/pyclone_output/{sample}.log"
    run:
        pur = samples.loc[wildcards.sample,'purity']
        dir = samples.loc[wildcards.sample, 'sample']
        outdir="data/pyclone_output/"+dir
        if not os.path.exists(outdir):
            os.makedirs(outdir)
        #import pdb; pdb.set_trace()
        #os.mkdir(dir,output[0])
        shell("PyClone run_analysis_pipeline --in_files {input} --working_dir {outdir} --tumour_contents {pur}  --density pyclone_binomial --min_cluster_size 2 > {output}")

rule clone_evol:
    input:
        "data/pyclone_output/{sample}/tables/loci.tsv"
    output:
        out_loci = "data/pyclone_output/{sample}/tables/loci_processed.tsv",
        out_tree = "data/clone_trees/{sample}_tree.txt"
    run:
        shell("Rscript scripts/clone_tree_inference_clonevol.R {input} {output.out_loci} {output.out_tree}")

rule get_clone_sequences:
    input:
        snv = "data/snv/{sample}_SNV.tsv.gz",
        loci = "data/pyclone_output/{sample}/tables/loci_processed.tsv",
        tree = "data/trees/{sample}_tree.txt"
    output:
        seqs = "data/codon_seqs/{sample}_clone_seqs.fas"
    shell:
        "Rscript scripts/get_clone_sequences.R {input.snv} {input.loci} {input.tree} {output.seqs}"

rule curate_clone_seqs:
    input:
        "data/seqs/{sample}_clone_seqs.fas"
    output:
        "data/seqs/{sample}_clone_seqs_curated.fas"
    shell:
        "sed 's/[\.]//g' {input} > {output}"

#rule convert_fasta_2_phylip:
#    input:
#        "data/codon_seqs/{sample}_clone_seqs_curated.fas"
#    output:
#        "data/nuc_seqs/{sample}_clone_seqs_curated.phy"
#    shell:
#        "python scripts/fasta2phylip.py -i {input} -o {output}"

#python scripts/ascbias.py -p data/seqs/{samples}_pyclone_seqs.phy -o data/seqs/{sample}_pyclone_seqs

rule build_tree:
    input: "data/seqs/{sample}_clone_seqs_curated.fas"
    output: "data/raxml_trees/{sample}_raxml.log"
    run:
        dir = samples.loc[wildcards.sample, 'sample']
        outdir= "data/raxml_trees/"+dir
        shell("raxml-ng --msa {input} --model GTR+G --force model_lh_impr --precision 12 --blmin 0.000000001 --prefix {outdir} > {output}")
