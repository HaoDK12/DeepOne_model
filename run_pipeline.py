# -*- coding: utf-8 -*-
"""
CRISPR Joint Prediction Pipeline: DeepOne + inDelphi
Description: This script dynamically generates and executes the dual-environment 
             Prediction pipeline of 1-bp insertion frequencies by bridging Python 3 (DeepOne) and Python 2.7 (inDelphi).
"""

import os
import sys
import subprocess
try:
    import configparser
except ImportError:
    import ConfigParser as configparser  # Compatibility fallback

# =====================================================================
# 1. LOAD CONFIGURATION FROM FILE
# =====================================================================

CONFIG_FILE = "config.ini"
if not os.path.exists(CONFIG_FILE):
    print("[CRITICAL ERROR] Configuration file '%s' not found in current directory." % CONFIG_FILE)
    sys.exit(1)

config = configparser.ConfigParser()
config.read(CONFIG_FILE)

try:
    SAMPLE_GENOMIC_SEQUENCE = config.get("SEQUENCE_AND_CELL", "SAMPLE_GENOMIC_SEQUENCE").strip()
    CELL_LINE = config.get("SEQUENCE_AND_CELL", "CELL_LINE").strip()
    
    PYTHON3_EXEC = config.get("DEEPONE_ENV", "PYTHON3_EXEC").strip()
    DEEPONE_DIR  = config.get("DEEPONE_ENV", "DEEPONE_DIR").strip()
    
    PYTHON2_EXEC = config.get("INDELPHI_ENV", "PYTHON2_EXEC").strip()
    INDELPHI_DIR = config.get("INDELPHI_ENV", "INDELPHI_DIR").strip()
except Exception as e:
    print("[CRITICAL ERROR] Failed to parse parameters from 'config.ini'. Details: %s" % str(e))
    sys.exit(1)

# 2. PARAMETER VALIDATION
# =====================================================================

print("=" * 80)
print("[Bioinformatics Pipeline] Initializing validation procedures...")
print("=" * 80)

# Validate Cell Line Intersection
SUPPORTED_CELLS = ['mESCs', 'HEK293', 'K562']
if CELL_LINE not in SUPPORTED_CELLS:
    print("[CRITICAL ERROR] Supported cell lines are limited to the intersection: %s" % str(SUPPORTED_CELLS))
    print("Your input '%s' is incompatible with either DeepOne or inDelphi." % CELL_LINE)
    sys.exit(1)

# In inDelphi, the string format for mESC is 'mESC'. Ensure exact mapping.
indelphi_cell_mapping = {"mESCs": "mESC", "HEK293": "HEK293", "K562": "K562"}
INDELPHI_CELLTYPE = indelphi_cell_mapping[CELL_LINE]

# Validate Directory Existence
if not os.path.isdir(DEEPONE_DIR):
    print("[CRITICAL ERROR] DeepOne repository directory not found at: %s" % DEEPONE_DIR)
    sys.path.exit(1)
if not os.path.isdir(INDELPHI_DIR):
    print("[CRITICAL ERROR] inDelphi repository directory not found at: %s" % INDELPHI_DIR)
    sys.exit(1)

# Verify Executable Paths via minimal checks
try:
    p3_check = subprocess.check_output([PYTHON3_EXEC, "--version"], stderr=subprocess.STDOUT, text=True)
    print("[INFO] Python 3 Environment Detected: %s" % p3_check.strip())
except Exception as e:
    print("[CRITICAL ERROR] Cannot execute Python 3 using '%s'. Check path or alias." % PYTHON3_EXEC)
    sys.exit(1)

try:
    p2_check = subprocess.check_output([PYTHON2_EXEC, "--version"], stderr=subprocess.STDOUT)
    print("[INFO] Python 2 Environment Detected: %s" % p2_check.strip())
except Exception as e:
    print("[CRITICAL ERROR] Cannot execute Python 2.7 using '%s'. Check path." % PYTHON2_EXEC)
    sys.exit(1)

print("[SUCCESS] Parameter validation passed. Proceeding to generation stage.\n")

# Intermediate and Final filenames
INTERMEDIATE_CSV = os.path.abspath("./DeepOne_intermediates.csv")
FINAL_TSV = os.path.abspath("./Final_Joint_Predictions.tsv")

# =====================================================================
# 3. STEP 01: DYNAMIC GENERATION & EXECUTION OF DEEPONE (Python 3)
# =====================================================================
print("-" * 80)
print("[STAGE 1] Preparing DeepOne (Templated +1bp insertion prediction in Python 3)...")
print("-" * 80)

# Standardize path slashes
abs_deepone_dir = os.path.abspath(DEEPONE_DIR).replace('\\', '/')

step01_script_content = f"""import sys
import os
import pandas as pd
import numpy as np

DEEPONE_PATH = '{abs_deepone_dir}'
if DEEPONE_PATH not in sys.path:
    sys.path.append(DEEPONE_PATH)
os.chdir(DEEPONE_PATH)

from utils.Energy_cal import EnergyCalculator
from utils.DL_model import get_model

def rev_comp(seq: str) -> str:
    comp = {{'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C', 'N': 'N'}}
    return ''.join(comp.get(b, 'N') for b in reversed(seq.upper()))

def scan_ngg_pam_sites(genomic_seq: str) -> list:
    sites = []
    genomic_seq = genomic_seq.upper()
    
    for i in range(len(genomic_seq) - 2):
        if genomic_seq[i+1:i+3] == "GG":
            if i >= 24 and (i + 6) <= len(genomic_seq):
                sites.append({{
                    'strand': '+',
                    'id': f"gRNA_plus_{{len(sites)+1}}",
                    'target_30bp': genomic_seq[i-24:i+6],
                    'global_cutsite': i
                }})
                
    rc_seq = rev_comp(genomic_seq)
    for i in range(len(rc_seq) - 2):
        if rc_seq[i+1:i+3] == "GG":
            if i >= 24 and (i + 6) <= len(rc_seq):
                sites.append({{
                    'strand': '-',
                    'id': f"gRNA_minus_{{len(sites)+1}}",
                    'target_30bp': rc_seq[i-24:i+6],
                    'global_cutsite': i
                }})
    return sites

def one_hot_encode(seq: str) -> np.ndarray:
    mapping = {{'A': [1,0,0,0], 'T': [0,1,0,0], 'C': [0,0,1,0], 'G': [0,0,0,1]}}
    return np.array([mapping.get(b, [0,0,0,0]) for b in seq])

def run_deepone_stage(target_seq, cell_line="{CELL_LINE}"):
    candidate_sites = scan_ngg_pam_sites(target_seq)
    if not candidate_sites:
        print("[Warning] No candidate NGG sites matched length restrictions inside the sequence.")
        return
        
    print(f"[DeepOne Core] Loaded {{len(candidate_sites)}} candidate NGG sgRNAs. Commencing deep modeling...")
    energy_calculator = EnergyCalculator()
    model = get_model(max_len_en=30, max_dep=4)
    weight_dir = os.path.join(DEEPONE_PATH, "utils/bestmodel_DeepOne", cell_line)
    weight_paths = [os.path.join(weight_dir, f) for f in os.listdir(weight_dir) if f.endswith(".h5")]
    
    results = []
    for site in candidate_sites:
        target = site['target_30bp']
        guide = target[4:24] + target[24:27]
        energy_feat = energy_calculator.get_energy_features_for_guides({{site['id']: [guide, guide]}})
        onehot = one_hot_encode(target).reshape(1, 30, 4)
        rna_dna = np.array([[energy_feat[site['id']]["RNA_DNA_eng"]]])
        
        preds = []
        for w in weight_paths:
            model.load_weights(w)
            preds.append(model.predict([onehot, rna_dna], verbose=0))
            
        results.append({{
            "ID": site['id'],
            "Target_30bp": target,
            "Target_gRNA": target[4:24],
            "Strand": site['strand'],
            "PAM": target[24:27],
            "GC%": round((target[4:24].count('G') + target[4:24].count('C')) / 20.0 * 100, 1),
            "DeepOne_score": round(float(np.mean(preds)), 2),
            "global_cutsite": site['global_cutsite']
        }})
        
    df = pd.DataFrame(results)
    output_path = r'{INTERMEDIATE_CSV}'
    df.to_csv(output_path, index=False)
    print(f"[DeepOne Core] Complete. Records exported successfully.")

if __name__ == '__main__':
    run_deepone_stage("{SAMPLE_GENOMIC_SEQUENCE}")
"""

# Write temporary Step 1 runner
step01_file = "./_generated_DeepOne_step01.py"
with open(step01_file, "w", encoding="utf-8") as f:
    f.write(step01_script_content)

print("[INFO] Executing DeepOne step01 sub-process...")
try:
    p1_run = subprocess.run([PYTHON3_EXEC, step01_file], check=True, capture_output=True, text=True)
    print(p1_run.stdout)
except subprocess.CalledProcessError as e:
    print("[CRITICAL ERROR] DeepOne execution failed!")
    print("--- Python 3 Error Output ---")
    print(e.stderr)
    sys.exit(1)

if not os.path.exists(INTERMEDIATE_CSV):
    print("[CRITICAL ERROR] Step 1 finished but intermediate CSV was not created. Check file writing permissions.")
    sys.exit(1)


# =====================================================================
# 4. STEP 02: DYNAMIC GENERATION & EXECUTION OF INDELPHI (Python 2.7)
# =====================================================================
print("-" * 80)
print("[STAGE 2] Preparing inDelphi (Indel spectrum profiling in Python 2.7)...")
print("-" * 80)

abs_indelphi_dir = os.path.abspath(INDELPHI_DIR).replace('\\', '/')

step02_script_content = f"""# -*- coding: utf-8 -*-
import sys
import os
import pandas as pd

INDELPHI_PATH = '{abs_indelphi_dir}'
if INDELPHI_PATH not in sys.path:
    sys.path.append(INDELPHI_PATH)
GLOBAL_GENOMIC_SEQ = '{SAMPLE_GENOMIC_SEQUENCE}'

import inDelphi

def main():
    input_csv = r'{INTERMEDIATE_CSV}'
    output_tsv = r'{FINAL_TSV}'
    
    if not os.path.exists(input_csv):
        print "[Error] Intermediate data file was not found."
        sys.exit(1)
        
    df = pd.read_csv(input_csv)
    print "[inDelphi Core] Initializing specific pre-trained weights for: {INDELPHI_CELLTYPE}..."
    inDelphi.init_model(celltype='{INDELPHI_CELLTYPE}')
    
    final_results = []
    for index, row in df.iterrows():
        g_cutsite = int(row['global_cutsite'])
        strand = row['Strand']
    
        if strand == '+':
            start_idx = max(0, g_cutsite - 37)
            end_idx = min(len(GLOBAL_GENOMIC_SEQ), g_cutsite + 30)
            context_seq = GLOBAL_GENOMIC_SEQ[start_idx:end_idx]
            cutsite_in_context = g_cutsite - start_idx
        else:
            def rev_comp(s):
                comp = {{'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C', 'N': 'N'}}
                return ''.join(comp.get(b, 'N') for b in reversed(s.upper()))
            rc_full = rev_comp(GLOBAL_GENOMIC_SEQ)
            rc_cutsite_pos = g_cutsite
            start_idx = max(0, rc_cutsite_pos - 37)
            end_idx = min(len(rc_full), rc_cutsite_pos + 30)
            context_seq = rc_full[start_idx:end_idx]
            cutsite_in_context = rc_cutsite_pos - start_idx
            
        try:
            pred_df, stats = inDelphi.predict(str(context_seq), int(cutsite_in_context))
            ins_1bp_freq = stats.get('1-bp ins frequency', 0.0)
            frameshift_freq = stats.get('Frameshift frequency', 1.0)
            
            ratio = ins_1bp_freq / frameshift_freq if frameshift_freq > 0 else 0.0
            
            final_results.append({{
                "ID": row['ID'],
                "Target_gRNA": row['Target_gRNA'],
                "Strand": row['Strand'],
                "PAM": row['PAM'],
                "GC%": row['GC%'],
                "DeepOne_score": row['DeepOne_score'],
                "inDelphi_1bp_Ins%": round(ins_1bp_freq, 2),
                "inDelphi_Frameshift%": round(frameshift_freq, 2),
                "Ins1bp_vs_Frameshift_Ratio": round(ratio, 4)
            }})
            print "Successfully computed landscape profile for site: " + str(row['ID'])
        except Exception as e:
            print "[Warning] Skipping site " + str(row['ID']) + " due to internal algorithm error: " + str(e)
            
    if not final_results:
        print "[Error] Matrix convergence failed. No rows output."
        sys.exit(1)
        
    out_df = pd.DataFrame(final_results)
    cols = ["ID", "Target_gRNA", "Strand", "PAM", "GC%", "DeepOne_score", "inDelphi_1bp_Ins%", "inDelphi_Frameshift%", "Ins1bp_vs_Frameshift_Ratio"]
    out_df = out_df[cols]
    out_df = out_df.sort_values(by="Ins1bp_vs_Frameshift_Ratio", ascending=False).reset_index(drop=True)
    out_df.to_csv(output_tsv, sep='\\t', index=False)
    print "[inDelphi Core] Batch profiling complete. Final ranked index generated."

if __name__ == '__main__':
    main()
"""

# Write temporary Step 2 runner
step02_file = "./_generated_inDelphi_step02.py"
with open(step02_file, "w", encoding="utf-8") as f:
    f.write(step02_script_content)

print("[INFO] Executing inDelphi step02 sub-process...")
try:
    # Run via Python 2 interpreter
    p2_run = subprocess.run([PYTHON2_EXEC, step02_file], check=True, capture_output=True, text=True)
    print(p2_run.stdout)
except subprocess.CalledProcessError as e:
    print("[CRITICAL ERROR] inDelphi execution failed!")
    print("--- Python 2.7 Error Output ---")
    print(e.stderr)
    sys.exit(1)

# =====================================================================
# 5. POST-PROCESSING & CLEANUP
# =====================================================================
print("=" * 80)
print("[PIPELINE COMPLETE SUCCESS]")
print("=" * 80)
print(">>> Intermediate File: %s" % INTERMEDIATE_CSV)
print(">>> FINAL RECOMMENDED TARGETS: %s" % FINAL_TSV)

# Automatically remove auto-generated script runners to maintain a clean workspace
for temp_file in [step01_file, step02_file]:
    if os.path.exists(temp_file):
        os.remove(temp_file)
