import sys
import pickle
import subprocess
from typing import Dict, List, Tuple

ENERGY_PARAMS_PATH = "/home/au735018/HTTPserver_Cas9design/Deepone/utils/energy_dics.pkl"

class EnergyCalculator:
    def __init__(self, energy_params_path=ENERGY_PARAMS_PATH):
        """初始化能量计算器"""
        self.RI_REV_NT_MAP = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C', 'N': 'N'}
        self.RI_MATCH_noGU = {
            'A': {'A': False, 'C': False, 'G': False, 'T': True},
            'C': {'A': False, 'C': False, 'G': True, 'T': False},
            'G': {'A': False, 'C': True, 'G': False, 'T': False},
            'T': {'A': True, 'C': False, 'G': False, 'T': False}
        }
        self.RI_DNA_DNA_NN = {'AA': {'TT': -1.00}, 'TT': {'AA': -1.00}, 'AT': {'TA': -0.88}, 'TA': {'AT': -0.58},
                 'CA': {'GT': -1.45}, 'TG': {'AC': -1.45}, 'GT': {'CA': -1.44}, 'AC': {'TG': -1.44},
                 'CT': {'GA': -1.28}, 'AG': {'TC': -1.28}, 'GA': {'CT': -1.30}, 'TC': {'AG': -1.30},
                 'CG': {'GC': -2.17}, 'GC': {'CG': -2.24}, 'GG': {'CC': -1.84}, 'CC': {'GG': -1.84}}


        ########## READ THE 2mer 3mer 4mer energies ###############
        self.RNA_DNA_internal_loop = {3: 3.2, 4: 3.555, 5: 3.725,
                         6: 3.975, 7: 4.16, 8: 4.33, 9: 4.495, 10: 4.6, 11: 4.7}
        self.POS_WGH = [1.80067099242007, 1.95666668400006, 1.90472004401173, 2.13047270152512, 1.37853848098249, 1.46460783730408, 1.0, 1.387220146823, 1.51401000729362, 1.98058344620751,
           1.87939168587699, 1.7222593588838, 2.02228445489326, 1.92692086621503, 2.08041972716723, 1.94496755678903, 2.14539112893591, 2.04277109036766, 2.24911493451185, 2.25]


        self.DNA_POS_WGH = [1.22245576981774, 1.24561578622024, 1.37883177517399, 1.39146340276523, 1.24308180746857, 1.09598194424544, 1.0, 1.11695025382169, 1.11589045394936, 1.22243614188218,
               1.21317477033274, 1.07125942316357, 1.25205871414019, 1.21445408158483, 1.20971491326295, 1.21076785001579, 1.2480898972246, 1.40301355270318, 1.41221084925493, 1.4]

        ######### pam correction parameters for pam-updated energy ##########
        self.pam_ratios = {"GGG": 1.0, "AGG": 1.0, "CGG": 1.0, "TGG": 1.0, "GAG": 0.9, "AAG": 0.9,
              "CAG": 0.9, "TAG": 0.9, "GGA": 0.8, "AGA": 0.8, "CGA": 0.8, "TGA": 0.8, 
              #"AGC":0.8,"AGA":1.0,"AGT":1.2,"AGG":1.0,##NG
              #"CGC":0.9,"CGA":1.0,"CGT":1.1,"CGG":1.0,
              #"GGC":0.9,"GGA":1.1,"GGT":1.2,"GGG":1.1,
              #"TGC":0.8,"TGA":0.9,"TGT":1.1,"TGG":1.0,
              "OTHERS": 0.0}
        self.REV_NT_MAP = {'-': '', 'a': 'T', 'A': 'T', 'c': 'G', 'C': 'G', 'g': 'C', 'G': 'C',
                't': 'A', 'T': 'A', 'u': 'A', 'U': 'A', 'n': 'N', 'N': 'N'}
        
        # 加载能量参数
        with open(energy_params_path, 'rb') as f:
            self.RNA_DNA = pickle.load(f)

        self.pam_ratio_count = 3
        

    def calcRNADNAenergy(self,guideSeq, otSeq, GU_allowed=False):
        guideSeq = guideSeq.upper()[:-3]
        seq = ''.join([self.RI_REV_NT_MAP[c] for c in otSeq[:-3]])

        spos = -1
        epos = -1

        energy = [0.0]*len(guideSeq)

        MATCH = self.RI_MATCH_noGU
        # if not GU_allowed:
        #    MATCH = RI_MATCH_noGU

        for i in range(len(seq)):
            if MATCH[guideSeq[i]][seq[i]]:
                if spos == -1:
                    spos = i
                epos = i

        i = spos
        while i < epos:
            j = i+1
            while MATCH[seq[j]][guideSeq[j]] == False:
                j = j+1

                if j > epos:
                    break
                
            loop_size = (j-i)-1
            eng_con = 0
            if loop_size < 3:
                # print i, j, loop_size, guideSeq[i:j+1], seq[i:j+1]
                eng_con = self.RNA_DNA[loop_size][guideSeq[i:j+1]][seq[i:j+1]]
                # if there is a stack in the beginning or end AU GU penalty is still needed
                if loop_size == 0:
                    if (i == spos and (guideSeq[i] == "T" or seq[i] == "T")) or (j == epos and (guideSeq[j] == "T" or seq[j] == "T")):
                        eng_con += 0.25
            elif loop_size <= 11:
                eng_con = float(self.RNA_DNA_internal_loop[loop_size]) + float(
                    self.RNA_DNA[0][guideSeq[i:i+2]][seq[i:i+2]]) + float(self.RNA_DNA[0][guideSeq[j-1:j+1]][seq[j-1:j+1]])

            for k in range(loop_size+1):
                energy[i+k] += eng_con/(loop_size+1)

            i = j

        return energy
    def calcDNAopening(self,otSeq):
        seq = otSeq.upper()[:-3]
        energy = [0.0]*len(seq)
        for i in range(1, len(seq)):
            energy[i] = float(self.RI_DNA_DNA_NN[seq[i-1]+seq[i]]
                            [self.RI_REV_NT_MAP[seq[i-1]]+self.RI_REV_NT_MAP[seq[i]]])
        return energy
    
    def get_rnafold_eng(self, seq_list):
        cmd = "RNAfold --noPS"
        input_data = "\n".join(seq_list)
        p = subprocess.Popen(cmd, shell=True, stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        stdout, _ = p.communicate(input=input_data.encode())
        results = stdout.decode().split("\n")
        energies = []
        for line in results:
          if "(" in line and ")" in line:
            try:
                eng = float(line.strip().split()[-1].replace("(", "").replace(")", ""))
                energies.append(eng)
            except:
                energies.append(0.0)
        return sum(energies)
    
    def get_eng(self,grna_seq, off_seq, score_func, GU_allowed=False, pos_weight=False, pam_corr=False, grna_folding=False, dna_opening=False, dna_pos_wgh=False):
        scores = score_func(grna_seq, off_seq, GU_allowed)

        if pos_weight:
            for i in range(len(scores)):
                if i < 21:
                    scores[-(i+1)] = self.POS_WGH[-(i+1)] * scores[-(i+1)]

        off = sum(scores) + 0.0
        off = (-1.0) * off

        if grna_folding:
            off += self.get_rnafold_eng(grna_seq[:20])

        # if aa_interaction:
        #    off += AA_NT_ENG["Arg"][off_seq[22]] + AA_NT_ENG["Arg"][off_seq[21]] + AA_NT_ENG["Lys"][REV_NT_MAP[off_seq[21]]] + AA_NT_ENG["Arg"][off_seq[17].replace("T","U")]

        if dna_opening:
            dna_scores = self.calcDNAopening(off_seq)
            if dna_pos_wgh:
                for i in range(len(dna_scores)):
                    if i < 21:
                        dna_scores[-(i+1)] = self.DNA_POS_WGH[-(i+1)] * \
                            dna_scores[-(i+1)]
            off += sum(dna_scores)

        if pam_corr:
            if off_seq[-self.pam_ratio_count:] in self.pam_ratios.keys():
                off = off * self.pam_ratios[off_seq[-self.pam_ratio_count:]]
            else:
                off = off * self.pam_ratios["OTHERS"]
        
        return off


    def rev_comp_seq(self,seq):
        s = ''
        for c in seq:
            s = self.REV_NT_MAP[c] + s
        return s


    def comp_seq(self,seq):
        s = ''
        for c in seq:
            s = s + self.REV_NT_MAP[c]
        return s
    def get_energy_features_for_guides(self, guide_seqs: Dict[str, List[str]]):
        """计算能量特征"""
        feature_out = {}
        for guide, guide_seq in guide_seqs.items():
            feature_out[guide] = {
                "RNA_DNA_eng": self.get_eng(guide_seq[0], guide_seq[1], self.calcRNADNAenergy, GU_allowed=False, pos_weight=True, pam_corr=True, grna_folding=True, dna_opening=True, dna_pos_wgh=False)
            }
        return feature_out
