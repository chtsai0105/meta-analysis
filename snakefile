import pandas as pd

configfile: "config.yaml"
sample_df = pd.read_csv(config['Metadata'], keep_default_na=False, na_values=['_'], comment="#")
path=config['Path']

target_list = list()
# target_list.extend(["{dir}/{acc}.fastq.gz".format(dir=path['fastq'], acc=acc) for acc in sample_df.loc[sample_df['Layout'] == 'single', 'Accession']])
# target_list.extend(["{dir}/{acc}_1.fastq.gz".format(dir=path['fastq'], acc=acc) for acc in sample_df.loc[sample_df['Layout'] == 'paired', 'Accession']])
# target_list.extend(["{dir}/{acc}_2.fastq.gz".format(dir=path['fastq'], acc=acc) for acc in sample_df.loc[sample_df['Layout'] == 'paired', 'Accession']])
target_list.extend(["{dir}/{acc}.fasta".format(dir=path['fasta'], acc=acc) for acc in sample_df['Accession']])


wildcard_constraints:
    acc = "[a-zA-Z]{3}\d+"

rule all:
    input:
        target_list

rule parallel_fastq_dump_single:
    output:
        "data/fastq/{acc}.fastq.gz"
    params:
        tmp_dir = "data/temp",
        output_dir = "data/fastq",
        params = ""
    threads:
        4
    envmodules:
        "parallel-fastq-dump"
    shell:
        """
        parallel-fastq-dump -T {params.tmp_dir} -O {params.output_dir} --threads {threads} --gzip {params.params} --sra-id {wildcards.acc}
        """

use rule parallel_fastq_dump_single as parallel_fastq_dump_paired with:
    output:
        R1 = "data/fastq/{acc}_1.fastq.gz",
        R2 = "data/fastq/{acc}_2.fastq.gz"
    params:
        tmp_dir = "data/temp",
        output_dir = "data/fastq",
        params = "--split-files"

rule bbduk_single:
    input:
        rules.parallel_fastq_dump_single.output
    output:
        "data/trim_fastq/{acc}.fastq.gz"
    params:
        "maq=10 qtrim=rl trimq=6 mlf=0.5 minlen=50"
    envmodules:
        "BBMap"
    shell:
        """
        bbduk.sh -Xmx1g in={input} out={output} {params}
        """

rule bbduk_paired:
    input:
        R1 = rules.parallel_fastq_dump_paired.output.R1,
        R2 = rules.parallel_fastq_dump_paired.output.R2
    output:
        R1 = "data/trim_fastq/{acc}_1.fastq.gz",
        R2 = "data/trim_fastq/{acc}_2.fastq.gz"
    params:
        "maq=10 qtrim=rl trimq=6 mlf=0.5 minlen=50"
    shell:
        """
        bbduk.sh -Xmx1g in={input.R1} in2={input.R2} out={output.R1} out2={output.R2} {params}
        """

rule vsearch_mergepairs:
    input:
        R1 = lambda wildcards: "data/trim_fastq/{acc}_1.fastq.gz".format(acc=sample_df.loc[(sample_df['Accession'] == wildcards.acc) & (sample_df['Layout'] == 'paired'), 'Accession'].squeeze()),
        R2 = lambda wildcards: "data/trim_fastq/{acc}_2.fastq.gz".format(acc=sample_df.loc[(sample_df['Accession'] == wildcards.acc) & (sample_df['Layout'] == 'paired'), 'Accession'].squeeze())
    output:
        "data/fasta/{acc}.fasta"
    envmodules:
        "vsearch"
    shell:
        """
        vsearch --fastq_mergepairs {input.R1} --reverse {input.R2} --fastq_allowmergestagger --fastq_maxdiffpct 10 --fastq_minovlen 10 --fastaout {output}
        """

rule reformat_to_fasta:
    input:
        lambda wildcards: "data/trim_fastq/{acc}.fastq.gz".format(acc=sample_df.loc[(sample_df['Accession'] == wildcards.acc) & (sample_df['Layout'] == 'single'), 'Accession'].squeeze())
    output:
        "data/fasta/{acc}.fasta"
    envmodules:
        "BBMap"
    shell:
        """
        reformat.sh -Xmx2g in={input} out={output}
        """
