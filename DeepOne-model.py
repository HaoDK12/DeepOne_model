import argparse
import sys
import os
import pickle
import re
import numpy as np
import pandas as pd
import tensorflow as tf
from typing import Dict, List, Tuple
from utils.Energy_cal import EnergyCalculator
from utils.DL_model import get_model
from utils.DL_model2 import build_model

# --------------------- setting parameter ---------------------
MAX_LEN_EN = 30
MAX_DEP = 4

# --------------------- functions ---------------------
def rev_comp(seq: str) -> str:
    comp = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C', 'N': 'N'}
    return ''.join([comp.get(base, 'N') for base in reversed(seq)])

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
    return [(strand, seq_id, target) for strand, seq_id, target in sites if len(target) == 30]

def one_hot_encode(seq: str) -> np.ndarray:
    mapping = {'A': [1,0,0,0], 'T': [0,1,0,0], 'C': [0,0,1,0], 'G': [0,0,0,1]}
    return np.array([mapping.get(base, [0,0,0,0]) for base in seq]).flatten()

def calculate_gc_content(seq: str) -> float:
    return (seq.count('G') + seq.count('C')) / len(seq)

def load_prediction_model(cell_line2: str):
    model_dirs = {
        "HEK293": "./utils/bestmodel_DeepOne/HEK293",
        "CHO": "./utils/bestmodel_DeepOne/CHO",
        "HAP1": "./utils/bestmodel_DeepOne/HAP1",
        "iPSC": "./utils/bestmodel_DeepOne/iPSC",
        "K562": "./utils/bestmodel_DeepOne/K562",
        "mESCs": "./utils/bestmodel_DeepOne/mESCs",
        "RPE-1": "./utils/bestmodel_DeepOne/RPE-1",
    }
    if cell_line2 not in model_dirs:
        raise ValueError(f"Unsupported cell line: {cell_line2}. Supported: {list(model_dirs.keys())}")
    model_dir = model_dirs[cell_line2]
    h5_files = [f for f in os.listdir(model_dir) if f.endswith('.h5')]
    if not h5_files:
        raise FileNotFoundError(f"No .h5 files found in {model_dir}")
    h5_path = os.path.join(model_dir, h5_files[0])
    with tf.device('/CPU:0'):
        model = get_model(MAX_LEN_EN, MAX_DEP)
        model.load_weights(h5_path)
        model.trainable = False
    return model

def process_sequence(input_seq: str, cell_line2: str) -> List[Dict]:
    pam_sites = find_pam_sites(input_seq)
    energy_calculator = EnergyCalculator()
    model = load_prediction_model(cell_line2)
    results = []
    for strand, seq_id, target in pam_sites:
        target2 = target[4:24] + target[-6:-3]
        energy_features = energy_calculator.get_energy_features_for_guides({seq_id: [target2, target2]})
        onehot_feature = one_hot_encode(target).reshape((-1, 30, 4))
        rna_dna_energy = np.array([[energy_features[seq_id]["RNA_DNA_eng"]]])
        with tf.device('/CPU:0'):
            prediction1 = model.predict([onehot_feature, rna_dna_energy], verbose=0)[0][0]
        model2 = build_model(30, 4, "CG")
        predictions = []
        model_paths = [
            "utils/bestmodel_CRISPRon/model1/variables/variables",
            "utils/bestmodel_CRISPRon/model2/variables/variables",
            "utils/bestmodel_CRISPRon/model3/variables/variables",
            "utils/bestmodel_CRISPRon/model4/variables/variables",
            "utils/bestmodel_CRISPRon/model5/variables/variables",
            "utils/bestmodel_CRISPRon/model6/variables/variables"
        ]
        for mp in model_paths:
            with tf.device('/CPU:0'):
                model2.load_weights(mp).expect_partial()
                for layer in model2.layers:
                    layer.set_weights([w.astype('float32') for w in layer.get_weights()])
                model2.trainable = False
                pred = model2.predict([onehot_feature, rna_dna_energy], verbose=0)
                predictions.append(pred)
        avg_prediction = np.mean(predictions)
        GC = calculate_gc_content(target2[:20])
        results.append({
            'ID': seq_id,
            'Target': target2,
            'Strand': strand,
            'Deepone_score': round(float(np.expm1(prediction1)), 1),
            'CRISPRon_score': round(float(avg_prediction), 1),
            'GC%': round(float(GC * 100), 1),
            'PAM': target[-6:-3] if re.match(r"[ATCG]GG", target[-6:-3]) else target[-6:-4]
        })
    return results

def main():
    parser = argparse.ArgumentParser(description="Predict gRNA efficiency with DeepOne and CRISPRon models.")
    parser.add_argument("--input_seq", type=str, required=True, help="Input DNA sequence.")
    parser.add_argument("--cell_line", type=str, required=True, help="Cell line name (e.g., HEK293, CHO, etc.)")
    parser.add_argument("--out_path", type=str, default="output.tsv", help="Output file path.")
    parser.add_argument("--prefix", type=str, default="", help="Optional prefix for each ID.")
    
    args = parser.parse_args()

    try:
        results = process_sequence(args.input_seq, args.cell_line)
        with open(args.out_path, 'w') as f:
            header = "ID\tTarget\tStrand\tDeepone_score\tCRISPRon_score\tGC%\tPAM"
            f.write(header + "\n")
            for res in results:
                res_id = f"{args.prefix}{res['ID']}"
                f.write(f"{res_id}\t{res['Target']}\t{res['Strand']}\t{res['Deepone_score']:.1f}\t"
                        f"{res['CRISPRon_score']:.1f}\t{res['GC%']:.1f}\t{res['PAM']}\n")
        print(f"^o^ Great! Job finished and results successfully saved, Check it in the {args.out_path} file ! Bye!")
    except Exception as e:
        print(f"Error: {str(e)}", file=sys.stderr)
        sys.exit(1)

if __name__ == "__main__":
    main()
