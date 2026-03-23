def make_cmd(case, control, ref, EDNAFULL):
    cmd = f"""#!/bin/bash
#SBATCH --account DREAM_CRISPR
#SBATCH -c 1
#SBATCH --mem 30g
#SBATCH --time 24:00:00


## remove synthesis errors
julia scripts/remove_synthesis_error.jl --spCas9 raw_reads/{case} --WT raw_reads/{control} --reference {ref} --outPath remove_synthesis/{case} --start 104 --fin 140
python3 get_indel_profile/get_indel_profile.py -r {ref} -i remove_synthesis/{case} -o indel_profile/{case} -s 103 -f 140 -m {EDNAFULL}
python3 scripts/remove_sequencing_error.py -r {ref} -i indel_profile/{case} -o remove_sequencing_error/{case}
python3 scripts/get_indel_count.py --ref {ref} --sampleName {case}
"""
    return cmd


def get_all_data(samplefile):
    ref = "refs/ref_on_144.fa"
    EDNAFULL = "/home/panxiaoguang/DREAM_CRISPR/TRAP2024/TRAP_on_merge/get_indel_profile/EDNAFULL"
    case_pair = []
    with open(samplefile, "r") as f:
        for line in f:
            sample_name = line.strip("\n")
            for celltype in ["HEK", "U2OS"]:
                for rep in ["1", "2"]:
                    case_name = f"{celltype}_on_{sample_name}-{rep}"
                    control_name = f"{celltype}_on_WT-{rep}"
                    case_pair.append((case_name, control_name))
    return case_pair, ref, EDNAFULL


def main():
    samplefile = "samples2"
    case_pair, ref, EDNAFULL = get_all_data(samplefile)
    for i, case in enumerate(case_pair):
        cmd = make_cmd(case[0], case[1], ref, EDNAFULL)
        with open(f"tk.{i}.sh", "w") as f:
            f.write(cmd)


if __name__ == "__main__":
    main()
