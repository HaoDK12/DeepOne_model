import pandas as pd
from pyfaidx import Fasta
import os
import fire
def parse_indels(indel_str):
    positions_1bp = []
    positions = []
    if indel_str:
        entries = indel_str.split('|')
        for entry in entries:
            pos_type, _ = entry.split(':')
            pos, size = pos_type.split(
                'I') if 'I' in pos_type else pos_type.split('D')
            if size == "1":
                positions_1bp.append(int(pos))
                positions.append(int(pos))
            else:
                if 'D' in pos_type:
                    for pos in range(int(pos), int(pos)+int(size)):
                        positions.append(int(pos))
                else:
                    positions.append(int(pos))
    return positions, positions_1bp


def filter_reads(ins_str, del_str, lftS, lftE, rgtS, rgtE):
    """ 根据indel位置过滤reads """
    insertions, insertions_1bp = parse_indels(ins_str)
    deletions, deletions_1bp = parse_indels(del_str)
    # 合并插入和删除位置
    indel_positions = set(insertions + deletions)
    indel_positions_1bp = set(insertions_1bp + deletions_1bp)
    # 检查是否有indel在1-25或32-37范围内
    if any(pos in range(lftS-1,lftE+1) or pos in range(rgtS, rgtE+1) for pos in indel_positions_1bp):
        # 检查是否有indel在26-31范围内
        if any(pos in range(lftE+1, rgtS) for pos in indel_positions):
            return True  # 保留reads
        return False  # 删除reads
    return True

def check_wildtype(ins_str, del_str, cut_start, cut_end):
    ## check if indel postion will across from 16-18
    insertions, _ = parse_indels(ins_str)
    deletions, _ = parse_indels(del_str)
    # 合并插入和删除位置
    indel_positions = set(insertions + deletions)
    if any(pos in range(cut_start,cut_end+1) for pos in indel_positions):
        return False
    return True

def parse_indels_size(indel_str):
    Size = []
    if indel_str:
        entries = indel_str.split('|')
        for entry in entries:
            pos_type, _ = entry.split(':')
            pos, size = pos_type.split(
                'I') if 'I' in pos_type else pos_type.split('D')
            Size.append(int(size))
    return sum(Size)


def parse_mutations_size(mut_str):
    entries = []
    if mut_str:
        entries = mut_str.split('|')
    return len(entries)

def restore_sequences(ref, query):
    ref_restored = []
    query_restored = []
    
    for r, q in zip(ref, query):
        if r != '-': 
            ref_restored.append(r)
            query_restored.append(q)
    
    return ''.join(ref_restored), ''.join(query_restored)

def make_frequency_table(fs, lftS, lftE, rgtS, rgtE):
    total_counts = 0
    Aligned_Sequence = []
    Reference_Sequence = []
    Unedited = []
    n_deleted = []
    n_inserted = []
    n_mutated = []
    Reads = []
    with open(fs, "r") as f:
        for line in f:
            if line.startswith("reference"):
                continue
            else:
                Une = False
                ref, query, count, ins, dels, mis = line.strip(
                    "\n").split("\t")
                if filter_reads(ins, dels, lftS, lftE, rgtS, rgtE):
                    newref,newquery = restore_sequences(ref, query)
                    Reference_Sequence.append(newref)
                    Aligned_Sequence.append(newquery)
                    n_inserted.append(parse_indels_size(ins))
                    n_deleted.append(parse_indels_size(dels))
                    n_mutated.append(parse_mutations_size(mis))
                    total_counts += int(count)
                    if check_wildtype(ins, dels, 26, 28):
                        Une = True
                    Unedited.append(Une)
                    Reads.append(int(count))
    df = pd.DataFrame({"Aligned_Sequence": Aligned_Sequence, "Reference_Sequence": Reference_Sequence,
                      "Unedited": Unedited, "n_deleted": n_deleted, "n_inserted": n_inserted, "n_mutated": n_mutated, "Reads": Reads})
    df = df.assign(readP=df['Reads']/total_counts * 100)
    df.columns = ["Aligned_Sequence", "Reference_Sequence", "Unedited",
                  "n_deleted", "n_inserted", "n_mutated", "#Reads", "%Reads"]
    df = df.sort_values(by="#Reads", ascending=False, ignore_index=True)
    return df

def main(reference="None",inputDir="None",outDir="None"):
    ref_reader = Fasta(reference)
    if not os.path.exists(outDir):
        os.makedirs(outDir)
    for record in ref_reader:
        inFile = os.path.join(inputDir, record.name + ".indelProfile.tsv")
        df = make_frequency_table(inFile, 1, 23, 34, 37)
        outFile = os.path.join(outDir, record.name + ".indelProfile.tsv")
        df.to_csv(outFile, sep="\t", index=False)

if __name__ == '__main__':
    fire.Fire(main)
