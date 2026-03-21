from CRISPResso2Align import global_align, read_matrix
from collections import Counter
import numpy as np
from pyfaidx import Fasta
import os

def ReadIn(fs):
    ## read a file line by line then push into a list
    seq_lst = []
    with open(fs, 'r') as f:
        for line in f:
            seq_lst.append(line.strip())
    return seq_lst

def guina(seqlist):
    ## read a list and then make a Counter
    seq = Counter(seqlist)
    return seq

def find_deletions(ref, query):
    deletions = []
    ref_pos = 0  # To keep track of the position in the original reference sequence

    i = 0
    while i < len(ref):
        if ref[i] == '-':
            # Skip the insertion positions in the reference
            i += 1
            continue

        if query[i] == '-':
            # A deletion is found
            del_start = ref_pos
            del_len = 0
            del_bases = []

            # Continue to count the length of the deletion
            while i < len(ref) and query[i] == '-':
                del_bases.append(ref[i])
                del_len += 1
                i += 1
                ref_pos += 1  # Only increment ref_pos inside the deletion loop

            deletions.append((del_start, del_len, ''.join(del_bases)))
        else:
            i += 1
            ref_pos += 1  # Increment ref_pos when not in deletion

    return deletions

def find_insertions(ref, query):
    insertions = []
    ref_pos = 0  # To keep track of the position in the original reference sequence

    i = 0
    while i < len(ref):
        if ref[i] == '-':
            # An insertion is found
            ins_start = ref_pos - 1  # Position in the original reference sequence, adjusted to the right side
            ins_len = 0
            ins_bases = []

            # Continue to count the length of the insertion
            while i < len(ref) and ref[i] == '-':
                if query[i] != '-':
                    ins_bases.append(query[i])
                    ins_len += 1
                i += 1

            if ins_len > 0:
                insertions.append((ins_start, ins_len, ''.join(ins_bases)))
        else:
            i += 1
            ref_pos += 1  # Increment ref_pos when not in insertion

    return insertions

def find_mismatches(ref, query):
    mismatches = []
    ref_pos = 0  # To keep track of the position in the original reference sequence

    for i in range(len(ref)):
        if ref[i] == '-':
            # Skip the insertion positions in the reference
            continue

        if query[i] != '-' and ref[i] != query[i]:
            # A mismatch is found
            mismatch_start = ref_pos
            mismatch_raw = ref[i]
            mismatch_mut = query[i]
            mismatches.append((mismatch_start, mismatch_raw, mismatch_mut))

        ref_pos += 1

    return mismatches

def calling(uniq_reads, reference, matrix,outf):
    outf.write("reference\tquery\tcount\tinsertions\tdeletions\tmismatches\n")
    for query, count in uniq_reads.items():
        seq, ref, _ = global_align(query, reference, matrix=matrix, gap_incentive=np.array([0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0],dtype=int), gap_open=-15, gap_extend=-2, )
        insertions = find_insertions(ref, seq)
        deletions = find_deletions(ref, seq)
        mismatches = find_mismatches(ref, seq)
        insert_fmt = "|".join([f"{st+1}I{le}:{bs}" for st,le,bs in insertions])
        delet_fmt = "|".join([f"{st+1}D{le}:{bs}" for st,le,bs in deletions])
        mism_fmt = "|".join([f"{st+1}:{raw}2{mut}" for st,raw,mut in mismatches])
        outf.write(f"{ref}\t{seq}\t{count}\t{insert_fmt}\t{delet_fmt}\t{mism_fmt}\n")

def main(reference, inpath, outpath, start, fin, maxtix_file):
    ref = Fasta(reference)
    if not os.path.exists(outpath):
        os.mkdir(outpath)
    for record in ref:
        ## 检查“inpath/{record.name}.filter.out"是否存在
        if os.path.exists(f"{inpath}/{record.name}.filter.out"):
            with open(f"{outpath}/{record.name}.indelProfile.tsv", 'w') as outf:
                matrix = read_matrix(maxtix_file)
                total_reads = ReadIn(f"{inpath}/{record.name}.filter.out")
                total_reads_count = guina(total_reads)
                calling(total_reads_count, str(ref[record.name][start:fin]), matrix, outf)
        else:
            print(f"{inpath}/{record.name}.filter.out not found")

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("-r", "--reference", help="reference fasta file", required=True)
    parser.add_argument("-i", "--inpath", help="input path", required=True)
    parser.add_argument("-o", "--outpath", help="output path", required=True)
    parser.add_argument("-s", "--start", type=int, help="start position", required=True)
    parser.add_argument("-f", "--fin", type=int, help="end position", required=True)
    parser.add_argument("-m", "--maxtix", help="maxtix file", default="EDNAFULL")

    args = parser.parse_args()
    main(args.reference, args.inpath, args.outpath, args.start, args.fin, args.maxtix)


