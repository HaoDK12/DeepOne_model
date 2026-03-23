from pyfaidx import Fasta
import os
def get_indel_eff(fs):
    total_reads = 0
    indel_reads = 0
    with open(fs,"r") as f:
        for line in f:
            if line.startswith("Aligned_Sequence"):
                continue
            Aligned_Sequence,Reference_Sequence,Unedited,n_deleted,n_inserted,n_mutated,Reads,ReadsP = line.strip("\n").split("\t")
            total_reads += int(Reads)
            if Unedited == "False" and (n_deleted != "0" or n_inserted != "0"):
                indel_reads += int(Reads)
    return total_reads,indel_reads

def main(ref,sampleName):
    ref_reader = Fasta(ref)
    outF = f"indelCounts/{sampleName}.indel_efficiency.tsv"
    with open(outF,"w") as out:
        out.write("gRNAID\tTotal_reads\tIndel_reads\n")
        for record in ref_reader:
            gID = record.name
            inF = f"remove_sequencing_error/{sampleName}/{gID}.indelProfile.tsv"
            if not os.path.exists(inF):
                out.write(f"{gID}\t0\t0\n")
            else:
                total_reads,indel_reads = get_indel_eff(inF)
                out.write(f"{gID}\t{total_reads}\t{indel_reads}\n")

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description='Get indel efficiency of each gRNA')
    parser.add_argument('--ref', required=True, help='Reference fasta file')
    parser.add_argument('--sampleName', required=True, help='Sample name')
    args = parser.parse_args()
    main(args.ref,args.sampleName)
