from typing import List, Union
import numpy.typing as npt
import numpy as np


def base_count(fastafile: str) -> List[int]:
    # 課題 1-1
    num_A = 0; num_T = 0; num_G = 0; num_C = 0
    with open(fastafile) as f:
        header = f.readline()
        #　メモリに配慮して書くこと、readlines()だとメモリ食う
        for line in f:  # こうすると1行だけ取り込むらしい
            for base in line:
                if base == "A": num_A += 1
                elif base == "T": num_T += 1
                elif base == "G": num_G += 1
                elif base == "C": num_C += 1

    return [num_A, num_T, num_G, num_C] # A, T, G, C

def gen_rev_comp_seq(fastafile: str) -> str:
    # 課題 1-2
    # 逆相補鎖とは、塩基を相補的に置換し、3'-5'の向きを逆にしたもの
    with open(fastafile, mode="r") as f:
        header = f.readline()
        reverse_array = ""
        for line in f:
            for base in line:
                if base == "A": reverse_array += "T"
                elif base == "T":reverse_array += "A"
                elif base == "G": reverse_array += "C"
                elif base == "C": reverse_array += "G"

    return reverse_array[::-1]

#AttributeError: `np.float_` was removed in the NumPy 2.0 release. Use `np.float64` instead.
def calc_gc_content(fastafile: str, window: int=1000, step: int=300) -> Union[npt.NDArray[np.float64], List[float]]:
    # 課題 1-3
    # 値を出力するところまで。matplotlibを使う部分は別途実装してください。
    with open(fastafile, mode="r") as f:
        header = f.readline()
        sequence = "".join(line.strip() for line in f)
        GC_contents = []

    for i in range (0, len(sequence)-window+1 , step): # 末尾の(window-1)長の部分にはwindowの端が入らない
        GC_counter = 0
        sub_sequence = sequence[i:i+window]
        GC_counter += sub_sequence.count("G")
        GC_counter += sub_sequence.count("C")
        GC_contents.append(float(GC_counter/window)*100)

    return GC_contents

def search_motif(fastafile: str, motif: str) -> List[str]:
    # 課題 1-4
    with open(fastafile, mode="r") as f:
        header = f.readline()
        sequence = "".join(line.strip() for line in f)
        len_mtf = len(motif)
        # motifで配列を分割して長さを調べる→部分配列が重複する場所にあるとき使えない
        # スライディングウィンドウみたいに調べないとだめそう？

    fwd_motifs_locations = []
    index = sequence.find(motif)
    while index != -1:
        fwd_motifs_locations.append(f"F{index+1}")
        index = sequence.find(motif, index+1)

    rev_sequence = gen_rev_comp_seq(fastafile)
    rev_motifs_locations = []
    index = rev_sequence.find(motif)
    while index != -1:
        rev_motifs_locations.append(f"R{str(len(rev_sequence)-index)}")
        index = rev_sequence.find(motif, index+1)

    return fwd_motifs_locations + rev_motifs_locations

def translate(fastafile: str) -> List[str]:
    # 課題 1-5
    def convert_DNA_to_RNA(seq: str) -> str:
        return str.replace(seq, "T", "U")

    def translate_frame_into_aminoacid(frame:str)-> str:
        codon_dict = {
            "AUG":"M", # Met (start)
            "UAA":"_",  "UAG":"_", "UGA":"_", # stop
            "UUU":"F", "UUC":"F", # Phe
            "UUA":"L", "UUG":"L", "CUU":"L", "CUC":"L", "CUA":"L", "CUG":"L", # Leu
            "AUU":"I", "AUC":"I", "AUA":"I", # Ile
            "GUU":"V", "GUC":"V", "GUA":"V", "GUG":"V", # Val
            "UCU":"S", "UCC":"S", "UCA":"S", "UCG":"S", "AGU":"S", "AGC":"S", # Ser
            "CCU":"P", "CCC": "P", "CCA": "P", "CCG": "P", # Pro
            "ACU": "T", "ACC": "T", "ACA": "T", "ACG": "T", # Thr
            "GCU": "A", "GCC": "A", "GCA": "A", "GCG": "A", # Ala
            "UAU":"Y", "UAC":"Y", # Tyr
            "CAU":"H", "CAC":"H", # His
            "CAA":"Q", "CAG":"Q", # Gln
            "AAU":"N", "AAC":"N", # Asn
            "AAA":"K", "AAG":"K", # Lys
            "GAU":"D", "GAC":"D", # Asp
            "GAA":"E", "GAG":"E", # Glu
            "UGU":"C", "UGC":"C", # Cys
            "UGG":"W", # Trp
            "CGU": "R", "CGC": "R", "CGA": "R", "CGG": "R", "AGA": "R", "AGG": "R", # Arg
            "GGU":"G", "GGC": "G", "GGA": "G", "GGG": "G" # Gly
        }
        return codon_dict.get(frame)

    def read_and_translate(seq: str)->str:
        result = []
        for k in range(3):
            isTranslating = False
            amino_seq = ""
            for i in range (k, len(seq)-2 , 3):
                amino_acid = translate_frame_into_aminoacid(seq[i:i+3])
                if (amino_acid == "M" and isTranslating == False):
                    isTranslating = True
                if (isTranslating == True):
                    amino_seq += amino_acid
                if (amino_acid == "_" and isTranslating == True):
                    isTranslating = False
                    result.append(amino_seq)
                    amino_seq = ""
            if (amino_seq != ""):
                result.append(amino_seq)
        return result

    with open(fastafile, mode="r") as f:
        header = f.readline()
        sequence = "".join(line.strip() for line in f)
        fwd_rna = convert_DNA_to_RNA(sequence)
        rev_rna = convert_DNA_to_RNA(gen_rev_comp_seq(fastafile))

    return read_and_translate(fwd_rna) + read_and_translate(rev_rna)

if __name__ == "__main__":
    #filepath = "data/NT_113952.1.fasta"
    filepath = "data/ATGCCGT.fasta"
    #filepath = "data/NC_000012.fasta"
    # 課題 1-1
    print(base_count(filepath))
    # 課題 1-2
    print(gen_rev_comp_seq(filepath))
    # 課題 1-3
    print(calc_gc_content(filepath,window=2, step=1)) #"data/ATGCCGT.fasta"
    #print(calc_gc_content(filepath, window=10, step=7)) #""data/NC_000012.fasta"
    # 課題 1-4
    #print(search_motif(filepath, "ATG"))
    print(search_motif(filepath,"CG")) #"data/ATGCCGT.fasta"
    #print(search_motif(filepath, motif="CTAGCAT")) #"data/NC_000012.fasta"
    # 課題 1-5
    print(translate(filepath))
