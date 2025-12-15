import argparse
import sys
import os
import re
import numpy as np
import tensorflow as tf
from typing import Dict, List, Tuple

from utils.Energy_cal import EnergyCalculator
from utils.DL_model import get_model

# --------------------- 参数配置 ---------------------
MAX_LEN_EN = 30
MAX_DEP = 4

# --------------------- 工具函数 ---------------------
def rev_comp(seq: str) -> str:
    comp = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C', 'N': 'N'}
    return ''.join(comp.get(b, 'N') for b in reversed(seq))

def find_pam_sites(seq: str) -> List[Tuple[str, str, str]]:
    sites = []
    seq = seq.upper()
    rc_seq = rev_comp(seq)

    for i in range(len(seq) - 2):
        if seq[i+1:i+3] == "GG" or (seq[i+1] == "G" and seq[i+2] in "ATCG"):
            if i >= 24 and i + 6 <= len(seq):
                target = seq[i-24:i+6]
                sites.append(('+', f"p_{len(sites)+1}", target))

    for i in range(len(rc_seq) - 2):
        if rc_seq[i+1:i+3] == "GG" or (rc_seq[i+1] == "G" and rc_seq[i+2] in "ATCG"):
            if i >= 24 and i + 6 <= len(rc_seq):
                target = rc_seq[i-24:i+6]
                sites.append(('-', f"m_{len(sites)+1}", target))

    return [(s, i, t) for s, i, t in sites if len(t) == 30]

def one_hot_encode(seq: str) -> np.ndarray:
    mapping = {'A':[1,0,0,0],'T':[0,1,0,0],'C':[0,0,1,0],'G':[0,0,0,1]}
    return np.array([mapping.get(b,[0,0,0,0]) for b in seq])

def calculate_gc(seq: str) -> float:
    return (seq.count('G') + seq.count('C')) / len(seq)

# --------------------- 权重加载 ---------------------
def load_weight_paths(model_dir: str) -> List[str]:
    if not os.path.isdir(model_dir):
        raise FileNotFoundError(model_dir)
    return [
        os.path.join(model_dir, f)
        for f in os.listdir(model_dir)
        if f.endswith(".h5")
    ]

def get_cellline_weights(cell_line: str) -> List[str]:
    base = "./utils/bestmodel_DeepOne"
    mapping = {
        "HEK293":"HEK293","CHO":"CHO","HAP1":"HAP1",
        "iPSC":"iPSC","K562":"K562","mESCs":"mESCs","RPE-1":"RPE-1"
    }
    if cell_line not in mapping:
        raise ValueError(f"Unsupported cell line: {cell_line}")
    return load_weight_paths(os.path.join(base, mapping[cell_line]))

def get_variant_weights(variant: str) -> List[str]:
    base = "./utils/bestmodel_DeepOne/variants"
    if variant == "SpCas9-NG":
        return load_weight_paths(os.path.join(base, "ng"))
    if variant == "SpG":
        return load_weight_paths(os.path.join(base, "spg"))
    return []

# --------------------- 主预测逻辑 ---------------------
def process_sequence(
    input_seq: str,
    cell_line: str,
    variants: List[str]
) -> List[Dict]:

    pam_sites = find_pam_sites(input_seq)
    energy_calculator = EnergyCalculator()
    model = get_model(MAX_LEN_EN, MAX_DEP)

    cell_weights = get_cellline_weights(cell_line)
    variant_weights = {v: get_variant_weights(v) for v in variants}

    results = []

    for strand, seq_id, target in pam_sites:
        guide = target[4:24] + target[-6:-3]
        energy = energy_calculator.get_energy_features_for_guides(
            {seq_id: [guide, guide]}
        )
        onehot = one_hot_encode(target).reshape(1,30,4)
        rna_dna = np.array([[energy[seq_id]["RNA_DNA_eng"]]])

        # DeepOne 主预测
        preds = []
        for w in cell_weights:
            model.load_weights(w)
            preds.append(model.predict([onehot, rna_dna], verbose=0))
        deepone_score = float(np.mean(preds))

        row = {
            "ID": seq_id,
            "Target": guide,
            "Strand": strand,
            "DeepOne_score": round(deepone_score, 2),
            "GC%": round(calculate_gc(guide[:20]) * 100, 1),
            "PAM": target[-6:-3]
        }

        # Variant 预测（可选）
        for v in variants:
            vp = []
            for w in variant_weights[v]:
                model.load_weights(w)
                vp.append(model.predict([onehot, rna_dna], verbose=0))
            row[f"DeepOne-{v}_score"] = round(float(np.mean(vp)), 2)

        results.append(row)

    return results

# --------------------- CLI ---------------------
def main():
    parser = argparse.ArgumentParser("DeepOne gRNA efficiency prediction")
    parser.add_argument("--input_seq", required=True)
    parser.add_argument("--cell_line", required=True)
    parser.add_argument(
        "--cas9 variant",
        choices=["None","SpCas9-NG","SpG","Both"],
        default="None",
        help="Optional Cas variants"
    )
    parser.add_argument("--out", default="output.tsv")
    parser.add_argument("--prefix", default="")

    args = parser.parse_args()

    variants = []
    if args.variant == "SpCas9-NG":
        variants = ["SpCas9-NG"]
    elif args.variant == "SpG":
        variants = ["SpG"]
    elif args.variant == "Both":
        variants = ["SpCas9-NG","SpG"]

    results = process_sequence(args.input_seq, args.cell_line, variants)

    with open(args.out, "w") as f:
        header = results[0].keys()
        f.write("\t".join(header) + "\n")
        for r in results:
            r["ID"] = args.prefix + r["ID"]
            f.write("\t".join(map(str, r.values())) + "\n")

    print(f"Done. Results saved to {args.out}")

if __name__ == "__main__":
    main()

