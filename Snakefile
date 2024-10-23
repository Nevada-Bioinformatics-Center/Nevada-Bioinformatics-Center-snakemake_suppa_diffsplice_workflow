from snakemake.utils import validate, min_version
import os
configfile: "config.yaml"

INPUT_DIR=config["input_data"]
ANNOT=config["annot"]
ioetypes=config["ioetypes"]
ioelist = ioetypes.split(",")
print(ioelist)

CWD = os.getcwd()

##### target rules #####
rule all:
    input:
        "output/iso_tpm.txt",
        "output/DT_diffSplice_iso.filtered05.dpsi",
        expand("output/DT_diffSplice_iso_ioe_{ioe}.filtered05.dpsi", ioe=ioelist),
        #"output/DT_diffSplice_iso_ioe.filtered05.dpsi",
        "output/WW_diffSplice_iso.filtered05.dpsi",
        "output/MU_diffSplice_iso.filtered05.dpsi",
        "output/WT_diffSplice_iso.filtered05.dpsi",

rule suppa_makeisotpm:
    input:
        INPUT_DIR,
    output:
        "output/iso_tpm.txt"
    log: "log/suppa/iso_tpm.log"
    threads: 1
    resources: time_min=320, mem_mb=40000, cpus=1
    conda: "suppa.yml"
    shell: "python ~/git/SUPPA/multipleFieldSelection.py -i {input}/*/quant.sf -k 1 -f 4 -o {output} &> {log}" 

rule suppa_genEvents:
    input:
        annot=ANNOT,
    output:
        "output/suppa_gen_events_trans.ioi",
    params:
        "output/suppa_gen_events_trans",
    log: "log/suppa/gen_events.log"
    threads: 1
    resources: time_min=320, mem_mb=40000, cpus=1
    conda: "suppa.yml"
    shell: "python ~/git/SUPPA/suppa.py generateEvents -f ioi -i {input.annot} -o {params} &> {log}"

#rule suppa_genEvents_ioe:
#    input:
#        annot=ANNOT,
#    output:
#        "output/suppa_gen_events_trans_{ioe}.ioe",
#    params:
#        out="output/suppa_gen_events_trans",
#        ioe="{ioe}",
#    log: "log/suppa/gen_events_{ioe}.log"
#    threads: 1
#    resources: time_min=320, mem_mb=40000, cpus=1
#    conda: "suppa.yml"
#    shell: "python ~/git/SUPPA/suppa.py generateEvents -f ioe -i {input.annot} -o {params.out} -e {params.ioe} &> {log}"

rule suppa_genEvents_ioe:
    input:
        annot=ANNOT,
    output:
        "output/suppa_gen_events_trans_A3_strict.ioe",
        "output/suppa_gen_events_trans_A5_strict.ioe",
        "output/suppa_gen_events_trans_AF_strict.ioe",
        "output/suppa_gen_events_trans_AL_strict.ioe",
        "output/suppa_gen_events_trans_MX_strict.ioe",
        "output/suppa_gen_events_trans_RI_strict.ioe",
        "output/suppa_gen_events_trans_SE_strict.ioe",
    params:
        out="output/suppa_gen_events_trans",
    log: "log/suppa/gen_events_ioe.log"
    threads: 1
    resources: time_min=320, mem_mb=40000, cpus=1
    conda: "suppa.yml"
    shell: "python ~/git/SUPPA/suppa.py generateEvents -f ioe -i {input.annot} -o {params.out} -e SE SS MX RI FL &> {log}"

rule suppa_psiPerIsoform:
    input:
        isotpm="output/iso_tpm.txt",
        annot=ANNOT,
    output:
        "output/iso_isoform.psi",
    params:
        "output/iso",
    log: "log/suppa/psiPerIsoform.log"
    threads: 1
    resources: time_min=320, mem_mb=40000, cpus=1
    conda: "suppa.yml"
    shell: 
        "python ~/git/SUPPA/suppa.py psiPerIsoform -g {input.annot} -e {input.isotpm} -o {params} &> {log}"

##DT MU vs WT
rule suppa_splitfile_events_DT:
    input:
        psi="output/iso_isoform.psi",
    output:
        cond1="output/MU_DT_iso.psi",
        cond2="output/WT_DT_iso.psi",
    params:
        cond1="MU-DT.rep1,MU-WW.rep2,MU-WW.rep3",
        cond2="WT-DT.rep1,WT-WW.rep2,WT-WW.rep3",
    log: "log/suppa/splitfile_events_DT.log"
    threads: 1
    resources: time_min=320, mem_mb=20000, cpus=1
    conda: "suppa.yml"
    shell: "/data/gpfs/assoc/inbre/hansvg/home/miniconda3/envs/suppa/bin/Rscript /data/gpfs/home/hvasquezgross/git/SUPPA/scripts/split_file.R {input} {params.cond1} {params.cond2} {output.cond1} {output.cond2} -e &> {log}"

rule suppa_splitfile_tpm_DT:
    input:
        psi="output/iso_tpm.txt",
    output:
        cond1="output/MU_DT_iso.tpm",
        cond2="output/WT_DT_iso.tpm",
    params:
        cond1="MU-DT.rep1,MU-WW.rep2,MU-WW.rep3",
        cond2="WT-DT.rep1,WT-WW.rep2,WT-WW.rep3",
    log: "log/suppa/splitfile_tpm_DT.log"
    threads: 1
    resources: time_min=320, mem_mb=20000, cpus=1
    conda: "suppa.yml"
    shell: "/data/gpfs/home/hvasquezgross/git/SUPPA/scripts/split_file.R {input} {params.cond1} {params.cond2} {output.cond1} {output.cond2} -i &> {log}"

rule suppa_diffsplice_DT:
    input:
        ioi="output/suppa_gen_events_trans.ioi",
        cond1="output/MU_DT_iso.psi",
        cond2="output/WT_DT_iso.psi",
        cond1tpm="output/MU_DT_iso.tpm",
        cond2tpm="output/WT_DT_iso.tpm",
    output:
        "output/DT_diffSplice_iso.dpsi.temp.0",
    params:
        "output/DT_diffSplice_iso",
    log: "log/suppa/diffSplice_DT.log"
    threads: 1
    resources: time_min=320, mem_mb=50000, cpus=1
    conda: "suppa.yml"
    shell: "python /data/gpfs/home/hvasquezgross/git/SUPPA/suppa.py diffSplice -m empirical -gc -i {input.ioi} -p {input.cond1} {input.cond2} -e {input.cond1tpm} {input.cond2tpm} -o {params} &> {log}"

rule suppa_filter_DT:
    input:
        "output/DT_diffSplice_iso.dpsi.temp.0",
    output:
        "output/DT_diffSplice_iso.filtered05.dpsi",
    params:
    threads: 1
    resources: time_min=320, mem_mb=20000, cpus=1
    shell: "awk '$3 <= 0.05 {{print $0}}' {input} > {output}"

rule suppa_diffsplice_DT_ioe:
    input:
        #ioe="output/suppa_gen_events_trans_{ioe}.ioe",
        ioe="output/suppa_gen_events_trans_{ioe}_strict.ioe",
        cond1="output/MU_DT_iso.psi",
        cond2="output/WT_DT_iso.psi",
        cond1tpm="output/MU_DT_iso.tpm",
        cond2tpm="output/WT_DT_iso.tpm",
    output:
        "output/DT_diffSplice_iso_ioe_{ioe}.dpsi.temp.0",
    params:
        "output/DT_diffSplice_iso_ioe_{ioe}",
    log: "log/suppa/diffSplice_DT_ioe_{ioe}.log"
    threads: 1
    resources: time_min=320, mem_mb=50000, cpus=1
    conda: "suppa.yml"
    shell: "python /data/gpfs/home/hvasquezgross/git/SUPPA/suppa.py diffSplice -m empirical -gc -i {input.ioe} -p {input.cond1} {input.cond2} -e {input.cond1tpm} {input.cond2tpm} -o {params} &> {log}"

rule suppa_filter_DT_ioe:
    input:
        "output/DT_diffSplice_iso_ioe_{ioe}.dpsi.temp.0",
    output:
        "output/DT_diffSplice_iso_ioe_{ioe}.filtered05.dpsi",
    params:
    threads: 1
    resources: time_min=320, mem_mb=20000, cpus=1
    shell: "awk '$3 <= 0.05 {{print $0}}' {input} > {output}"

##WW MU vs WT
rule suppa_splitfile_events_WW:
    input:
        psi="output/iso_isoform.psi",
    output:
        cond1="output/MU_WW_iso.psi",
        cond2="output/WT_WW_iso.psi",
    params:
        cond1="MU-WW.rep1,MU-WW.rep2,MU-WW.rep3",
        cond2="WT-WW.rep1,WT-WW.rep2,WT-WW.rep3",
    log: "log/suppa/splitfile_events_WW.log"
    threads: 1
    resources: time_min=320, mem_mb=20000, cpus=1
    conda: "suppa.yml"
    shell: "/data/gpfs/assoc/inbre/hansvg/home/miniconda3/envs/suppa/bin/Rscript /data/gpfs/home/hvasquezgross/git/SUPPA/scripts/split_file.R {input} {params.cond1} {params.cond2} {output.cond1} {output.cond2} -e &> {log}"

rule suppa_splitfile_tpm_WW:
    input:
        psi="output/iso_tpm.txt",
    output:
        cond1="output/MU_WW_iso.tpm",
        cond2="output/WT_WW_iso.tpm",
    params:
        cond1="MU-WW.rep1,MU-WW.rep2,MU-WW.rep3",
        cond2="WT-WW.rep1,WT-WW.rep2,WT-WW.rep3",
    log: "log/suppa/splitfile_tpm_WW.log"
    threads: 1
    resources: time_min=320, mem_mb=20000, cpus=1
    conda: "suppa.yml"
    shell: "/data/gpfs/home/hvasquezgross/git/SUPPA/scripts/split_file.R {input} {params.cond1} {params.cond2} {output.cond1} {output.cond2} -i &> {log}"

rule suppa_diffsplice_WW:
    input:
        ioi="output/suppa_gen_events_trans.ioi",
        cond1="output/MU_WW_iso.psi",
        cond2="output/WT_WW_iso.psi",
        cond1tpm="output/MU_WW_iso.tpm",
        cond2tpm="output/WT_WW_iso.tpm",
    output:
        "output/WW_diffSplice_iso.dpsi.temp.0",
    params:
        "output/WW_diffSplice_iso",
    log: "log/suppa/diffSplice_WW.log"
    threads: 1
    resources: time_min=320, mem_mb=20000, cpus=1
    conda: "suppa.yml"
    shell: "python /data/gpfs/home/hvasquezgross/git/SUPPA/suppa.py diffSplice -m empirical -gc -i {input.ioi} -p {input.cond1} {input.cond2} -e {input.cond1tpm} {input.cond2tpm} -o {params} &> {log}"

rule suppa_filter_WW:
    input:
        "output/WW_diffSplice_iso.dpsi.temp.0",
    output:
        "output/WW_diffSplice_iso.filtered05.dpsi",
    params:
    threads: 1
    resources: time_min=320, mem_mb=20000, cpus=1
    shell: "awk '$3 <= 0.05 {{print $0}}' {input} > {output}"


##WT DT vs WW
rule suppa_diffsplice_WT:
    input:
        ioi="output/suppa_gen_events_trans.ioi",
        cond1="output/WT_WW_iso.psi",
        cond2="output/WT_DT_iso.psi",
        cond1tpm="output/WT_WW_iso.tpm",
        cond2tpm="output/WT_DT_iso.tpm",
    output:
        "output/WT_diffSplice_iso.dpsi.temp.0",
    params:
        "output/WT_diffSplice_iso",
    log: "log/suppa/diffSplice_WT.log"
    threads: 1
    resources: time_min=320, mem_mb=20000, cpus=1
    conda: "suppa.yml"
    shell: "python /data/gpfs/home/hvasquezgross/git/SUPPA/suppa.py diffSplice -m empirical -gc -i {input.ioi} -p {input.cond1} {input.cond2} -e {input.cond1tpm} {input.cond2tpm} -o {params} &> {log}"

rule suppa_filter_WT:
    input:
        "output/WT_diffSplice_iso.dpsi.temp.0",
    output:
        "output/WT_diffSplice_iso.filtered05.dpsi",
    params:
    threads: 1
    resources: time_min=320, mem_mb=20000, cpus=1
    shell: "awk '$3 <= 0.05 {{print $0}}' {input} > {output}"

##MU DT vs WW
rule suppa_diffsplice_MU:
    input:
        ioi="output/suppa_gen_events_trans.ioi",
        cond1="output/MU_WW_iso.psi",
        cond2="output/MU_DT_iso.psi",
        cond1tpm="output/MU_WW_iso.tpm",
        cond2tpm="output/MU_DT_iso.tpm",
    output:
        "output/MU_diffSplice_iso.dpsi.temp.0",
    params:
        "output/MU_diffSplice_iso",
    log: "log/suppa/diffSplice_MU.log"
    threads: 1
    resources: time_min=320, mem_mb=20000, cpus=1
    conda: "suppa.yml"
    shell: "python /data/gpfs/home/hvasquezgross/git/SUPPA/suppa.py diffSplice -m empirical -gc -i {input.ioi} -p {input.cond1} {input.cond2} -e {input.cond1tpm} {input.cond2tpm} -o {params} &> {log}"

rule suppa_filter_MU:
    input:
        "output/MU_diffSplice_iso.dpsi.temp.0",
    output:
        "output/MU_diffSplice_iso.filtered05.dpsi",
    params:
    threads: 1
    resources: time_min=320, mem_mb=20000, cpus=1
    shell: "awk '$3 <= 0.05 {{print $0}}' {input} > {output}"
