import pandas as pd

configfile: "config.yaml"
sample_df = pd.read_csv(config['Metadata'], keep_default_na=False, na_values=['_'], comment="#")
path=config['Path']

target_list = list()
# target_list.extend(["{dir}/{acc}.fastq.gz".format(dir=path['fastq'], acc=acc) for acc in sample_df.loc[sample_df['Layout'] == 'single', 'Accession']])
# target_list.extend(["{dir}/{acc}_1.fastq.gz".format(dir=path['fastq'], acc=acc) for acc in sample_df.loc[sample_df['Layout'] == 'paired', 'Accession']])
# target_list.extend(["{dir}/{acc}_2.fastq.gz".format(dir=path['fastq'], acc=acc) for acc in sample_df.loc[sample_df['Layout'] == 'paired', 'Accession']])
# target_list.extend(["{dir}/{kingdom}/{acc}.fasta".format(dir=path['fasta'], kingdom=v['Kingdom'], acc=v['Accession']) for _, v in sample_df[['Accession', 'Kingdom']].iterrows()])
target_list.extend(["{dir}/{kingdom}/{acc}.fasta".format(dir=path['rrna_prediction'], kingdom=v['Kingdom'], acc=v['Accession']) for _, v in sample_df[['Accession', 'Kingdom']].iterrows()])
# target_list.extend(["{dir}/{kingdom}.fasta".format(dir=path['rrna_prediction'], kingdom=kingdom) for kingdom in sample_df['Kingdom'].unique()])


wildcard_constraints:
    acc = "[a-zA-Z]{1,3}\d+",
    kingdom = "Bacteria|Fungi|Eukaryotes"

rule all:
    input:
        target_list

rule parallel_fastq_dump_single:
    output:
        temp("{dir}/{{acc}}.fastq.gz".format(dir=path['fastq']))
    params:
        tmp_dir = path['temp'],
        output_dir = path['fastq'],
        params = ""
    threads:
        8
    group:
        "single-end_preprossing"
    resources:
        time="1-00:00:00",
        mem_mb=lambda wildcards, input, attempt: 2500 * (1 + attempt)
    envmodules:
        "parallel-fastq-dump"
    shell:
        """
        parallel-fastq-dump -T {params.tmp_dir} -O {params.output_dir} --threads {threads} --gzip {params.params} --sra-id {wildcards.acc}
        """

use rule parallel_fastq_dump_single as parallel_fastq_dump_paired with:
    output:
        R1 = temp("{dir}/{{acc}}_1.fastq.gz".format(dir=path['fastq'])),
        R2 = temp("{dir}/{{acc}}_2.fastq.gz".format(dir=path['fastq']))
    params:
        tmp_dir = path['temp'],
        output_dir = path['fastq'],
        params = "--split-files"
    group:
        "paired-end_preprossing"

rule bbduk_single:
    input:
        rules.parallel_fastq_dump_single.output
    output:
        temp("{dir}/{{acc}}.fastq.gz".format(dir=path['trim_fastq']))
    params:
        "maq=10 qtrim=rl trimq=6 mlf=0.5 minlen=50"
    group:
        "single-end_preprossing"
    resources:
        time="1-00:00:00",
        mem_mb=lambda wildcards, input, attempt: min(max((input.size // 1000000) * (1 + attempt), 2000), 20000)
    envmodules:
        "BBMap"
    shell:
        """
        bbduk.sh -Xmx1g in={input} out={output} {params} && [[ -s {output} ]]
        """

rule bbduk_paired:
    input:
        R1 = rules.parallel_fastq_dump_paired.output.R1,
        R2 = rules.parallel_fastq_dump_paired.output.R2
    output:
        R1 = temp("{dir}/{{acc}}_1.fastq.gz".format(dir=path['trim_fastq'])),
        R2 = temp("{dir}/{{acc}}_2.fastq.gz".format(dir=path['trim_fastq']))
    params:
        "maq=10 qtrim=rl trimq=6 mlf=0.5 minlen=50"
    group:
        "paired-end_preprossing"
    resources:
        time="1-00:00:00",
        mem_mb=lambda wildcards, input, attempt: min(max((input.size // 1000000) * (1 + attempt), 2000), 20000)
    envmodules:
        "BBMap"
    shell:
        """
        bbduk.sh -Xmx1g in={input.R1} in2={input.R2} out={output.R1} out2={output.R2} {params} && [[ -s {output.R1} ]] && [[ -s {output.R2} ]]
        """

rule vsearch_mergepairs:
    input:
        R1 = lambda wildcards: "data/trim_fastq/{acc}_1.fastq.gz".format(acc=sample_df.loc[(sample_df['Accession'] == wildcards.acc) & (sample_df['Layout'] == 'paired'), 'Accession'].squeeze()),
        R2 = lambda wildcards: "data/trim_fastq/{acc}_2.fastq.gz".format(acc=sample_df.loc[(sample_df['Accession'] == wildcards.acc) & (sample_df['Layout'] == 'paired'), 'Accession'].squeeze())
    output:
        merged = temp("{dir}/{{acc}}_merged.fasta".format(dir=path['temp'])),
        notmerged_R1 = temp("{dir}/{{acc}}_notmerged_R1.fastq".format(dir=path['temp'])),
        notmerged_R2 = temp("{dir}/{{acc}}_notmerged_R2.fastq".format(dir=path['temp']))
    group:
        "paired-end_preprossing"
    resources:
        time="1-00:00:00",
        mem_mb=lambda wildcards, input, attempt: min(max((input.size // 1000000) * (1 + attempt), 2000), 20000)
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
        notmerged_interleaved = temp("{dir}/{{kingdom}}/{{acc}}_notmerged_interleaved.fasta".format(dir=path['fasta'])),
        intermediate = temp("{dir}/{{kingdom}}/{{acc}}_temp.fasta".format(dir=path['fasta'])),
        concat_fasta = temp("{dir}/{{kingdom}}/{{acc}}.fasta".format(dir=path['fasta']))
    group:
        "paired-end_preprossing"
    resources:
        time="1-00:00:00",
        mem_mb=lambda wildcards, input, attempt: min(max((input.size // 1000000) * (1 + attempt), 2000), 20000)
    envmodules:
        "BBMap"
    shell:
        """
        reformat.sh -Xmx2g in={input.notmerged_R1} in2={input.notmerged_R2} out={output.notmerged_interleaved}
        awk 'BEGIN{{OFS=FS=" "}}{{if(/^>/){{CUR=$1;{{if(CUR==PRE){{NUM++}}else{{NUM=1}}}};$1="";print CUR"."NUM $0;PRE=CUR}}else{{print $0}}}}' {output.notmerged_interleaved} > {output.intermediate}
        cat {input.merged} {output.intermediate} > {output.concat_fasta}
        """

rule reformat_to_fasta:
    input:
        lambda wildcards: "data/trim_fastq/{acc}.fastq.gz".format(acc=sample_df.loc[(sample_df['Accession'] == wildcards.acc) & (sample_df['Layout'] == 'single'), 'Accession'].squeeze())
    output:
        temp("{dir}/{{kingdom}}/{{acc}}.fasta".format(dir=path['fasta']))
    group:
        "single-end_preprossing"
    resources:
        time="1-00:00:00",
        mem_mb=lambda wildcards, input, attempt: min(max((input.size // 1000000) * (1 + attempt), 2000), 20000)
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
        marker_fasta = temp("{dir}/{{kingdom}}/temp/{{acc}}.fasta".format(dir=path['rrna_prediction']))
    params:
        lambda wildcards: sample_df.loc[sample_df['Accession'] == wildcards.acc, 'barnnap_key'].squeeze()
    threads:
        8
    resources:
        time="7-00:00:00",
        mem_mb=lambda wildcards, input, attempt: min(max((input.size // 1000000) * 5 * (1 + attempt * 1), 2000), 50000)
    conda:
        "envs/barrnap.yaml"
    shell:
        """
        barrnap --threads {threads} --kingdom {params} --reject 0.01 --lencutoff 0.01 --quiet --outseq {output.marker_fasta} {input} > /dev/null
        """

rule remove_duplicate:
    input: rules.barnapp.output.marker_fasta
    output: temp("{dir}/{{kingdom}}/temp/{{acc}}_nodup.fasta".format(dir=path['rrna_prediction']))
    group:
        "post-preprossing"
    conda:
        "envs/post_process.yaml"
    shell:
        """
        seqkit rmdup -n < {input} > {output}
        """

rule remove_N:
    input: rules.remove_duplicate.output
    output: "{dir}/{{kingdom}}/{{acc}}.fasta".format(dir=path['rrna_prediction'])
    group:
        "post-preprossing"
    conda:
        "envs/post_process.yaml"
    shell:
        """
        ./scripts/removeNfromfasta.py {input} > {output}
        """

rule combine_all_markers:
    input:
        lambda wildcards: expand("{dir}/{{kingdom}}/{acc}.fasta", dir=path['rrna_prediction'], acc=sample_df.loc[(sample_df['Kingdom'] == wildcards.kingdom), 'Accession'].squeeze())
    output:
        merged_fasta = "{dir}/{{kingdom}}.fasta".format(dir=path['rrna_prediction'])
    singularity:
        "envs/micca.simg"
    shell:
        """
        micca merge -i {input} -o {output.merged_fasta} -f fasta
        """
