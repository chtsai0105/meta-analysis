import pandas as pd

configfile: "config.yaml"
sample_df = pd.read_csv(config['Metadata'], keep_default_na=False, na_values=['_'], comment="#")
path=config['Path']

target_list = list()
# target_list.extend(["{dir}/{acc}.fastq.gz".format(dir=path['fastq'], acc=acc) for acc in sample_df.loc[sample_df['Layout'] == 'single', 'Accession']])
# target_list.extend(["{dir}/{acc}_1.fastq.gz".format(dir=path['fastq'], acc=acc) for acc in sample_df.loc[sample_df['Layout'] == 'paired', 'Accession']])
# target_list.extend(["{dir}/{acc}_2.fastq.gz".format(dir=path['fastq'], acc=acc) for acc in sample_df.loc[sample_df['Layout'] == 'paired', 'Accession']])
# target_list.extend(["{dir}/{acc}.fasta".format(dir=path['fasta'], acc=acc) for acc in sample_df['Accession']])
target_list.extend(["{dir}/{kingdom}/{acc}.fasta".format(dir=path['rrna_prediction'], kingdom=v['Kingdom'], acc=v['Accession']) for _, v in sample_df[['Accession', 'Kingdom']].iterrows()])
target_list.extend(["{dir}/{kingdom}.fasta".format(dir=path['rrna_prediction'], kingdom=kingdom) for kingdom in sample_df['Kingdom'].unique()])


wildcard_constraints:
    acc = "[a-zA-Z]{3}\d+",
    kingdom = "Bacteria|Fungi|Eukaryotes"

rule all:
    input:
        target_list

rule parallel_fastq_dump_single:
    output:
        "{dir}/{{acc}}.fastq.gz".format(dir=path['fastq'])
    params:
        tmp_dir = path['temp'],
        output_dir = path['fastq'],
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
        R1 = "{dir}/{{acc}}_1.fastq.gz".format(dir=path['fastq']),
        R2 = "{dir}/{{acc}}_2.fastq.gz".format(dir=path['fastq'])
    params:
        tmp_dir = path['temp'],
        output_dir = path['fastq'],
        params = "--split-files"

rule bbduk_single:
    input:
        rules.parallel_fastq_dump_single.output
    output:
        "{dir}/{{acc}}.fastq.gz".format(dir=path['trim_fastq'])
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
        R1 = "{dir}/{{acc}}_1.fastq.gz".format(dir=path['trim_fastq']),
        R2 = "{dir}/{{acc}}_2.fastq.gz".format(dir=path['trim_fastq'])
    params:
        "maq=10 qtrim=rl trimq=6 mlf=0.5 minlen=50"
    envmodules:
        "BBMap"
    shell:
        """
        bbduk.sh -Xmx1g in={input.R1} in2={input.R2} out={output.R1} out2={output.R2} {params}
        """

rule vsearch_mergepairs:
    input:
        R1 = lambda wildcards: "data/trim_fastq/{acc}_1.fastq.gz".format(acc=sample_df.loc[(sample_df['Accession'] == wildcards.acc) & (sample_df['Layout'] == 'paired'), 'Accession'].squeeze()),
        R2 = lambda wildcards: "data/trim_fastq/{acc}_2.fastq.gz".format(acc=sample_df.loc[(sample_df['Accession'] == wildcards.acc) & (sample_df['Layout'] == 'paired'), 'Accession'].squeeze())
    output:
        merged = temp("{dir}/{{acc}}_merged.fasta".format(dir=path['temp'])),
        notmerged_R1 = temp("{dir}/{{acc}}_notmerged_R1.fastq".format(dir=path['temp'])),
        notmerged_R2 = temp("{dir}/{{acc}}_notmerged_R2.fastq".format(dir=path['temp']))
    envmodules:
        "vsearch"
    shell:
        """
        vsearch --fastq_mergepairs {input.R1} --reverse {input.R2} --fastq_allowmergestagger --fastq_maxdiffpct 10 --fastq_minovlen 10 \
        --fastaout {output.merged} --fastqout_notmerged_fwd {output.notmerged_R1} --fastqout_notmerged_rev {output.notmerged_R2}
        """

rule concat_merged_and_unmerged_pairs:
    input:
        merged = rules.vsearch_mergepairs.output.merged,
        notmerged_R1 = rules.vsearch_mergepairs.output.notmerged_R1,
        notmerged_R2 = rules.vsearch_mergepairs.output.notmerged_R2,
    output:
        "{dir}/{{kingdom}}/{{acc}}.fasta".format(dir=path['fasta'])
    envmodules:
        "BBMap"
    shell:
        """
        reformat.sh -Xmx2g in={input.notmerged_R1} in2={input.notmerged_R2} out={wildcards.acc}_notmerged_interleaved.fasta
        awk 'BEGIN{{OFS=FS=" "}}{{if(/^>/){{CUR=$1;{{if(CUR==PRE){{NUM++}}else{{NUM=1}}}};$1="";print CUR"."NUM $0;PRE=CUR}}else{{print $0}}}}' {wildcards.acc}_notmerged_interleaved.fasta > {wildcards.acc}_temp.fasta
        cat {input.merged} {wildcards.acc}_temp.fasta > {output}
        rm -f {wildcards.acc}_notmerged_interleaved.fasta {wildcards.acc}_temp.fasta
        """

rule reformat_to_fasta:
    input:
        lambda wildcards: "data/trim_fastq/{acc}.fastq.gz".format(acc=sample_df.loc[(sample_df['Accession'] == wildcards.acc) & (sample_df['Layout'] == 'single'), 'Accession'].squeeze())
    output:
        "{dir}/{{kingdom}}/{{acc}}.fasta".format(dir=path['fasta'])
    envmodules:
        "BBMap"
    shell:
        """
        reformat.sh -Xmx2g in={input} out={output}
        """

rule barnapp:
    input:
        "{dir}/{{kingdom}}/{{acc}}.fasta".format(dir=path['fasta'])
    output:
        fai = temp("{dir}/{{kingdom}}/{{acc}}.fasta.fai".format(dir=path['fasta'])),
        marker_fasta = "{dir}/{{kingdom}}/{{acc}}.fasta".format(dir=path['rrna_prediction'])
    params:
        lambda wildcards: sample_df.loc[sample_df['Accession'] == wildcards.acc, 'barnnap_key'].squeeze()
    conda:
        "envs/barrnap.yaml"
    shell:
        """
        barrnap --kingdom {params} --reject 0.01 --lencutoff 0.01 --outseq {output.marker_fasta} {input}
        """

rule combine_all_markers:
    input:
        expand("{dir}/{{kingdom}}/{acc}.fasta", dir=path['rrna_prediction'], acc=sample_df['Accession'])
    output:
        "{dir}/{{kingdom}}.fasta".format(dir=path['rrna_prediction'])
    shell:
        """
        cat {input} > {output}
        """