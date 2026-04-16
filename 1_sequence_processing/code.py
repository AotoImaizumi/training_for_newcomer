from typing import List, Union
import numpy.typing as npt
import numpy as np

def base_count(fastafile: str) -> List[int]:
    # 課題 1-1
    num_A = 0; num_T = 0; num_G = 0; num_C = 0
    with open(fastafile) as f:
        header = f.readline()
        #　メモリに配慮して書くこと
        for line in f:
            for base in line:
                if base == "A": num_A += 1
                elif base == "T": num_T += 1
                elif base == "G": num_G += 1
                elif base == "C": num_C += 1

    return [num_A, num_T, num_G, num_C] # A, T, G, C

def gen_rev_comp_seq(fastafile: str) -> str:
    # 課題 1-2
    with open(fastafile, mode="r") as f:
        header = f.readline()
        reverse_array = ""
        for line in f:
            for base in line:
                if base == "A": reverse_array += "T"
                elif base == "C":reverse_array += "A"
                elif base == "T": reverse_array += "G"
                elif base == "G": reverse_array += "C"


    return reverse_array

def calc_gc_content(fastafile: str, window: int=1000, step: int=300) -> Union[npt.NDArray[np.float_], List[float]]:
    # 課題 1-3
    # 値を出力するところまで。matplotlibを使う部分は別途実装してください。
    with open(fastafile, mode="r") as f:
        header = f.readline()
        array = f.readlins()
    gc_content = []
    return []

def search_motif(fastafile: str, motif: str) -> List[str]:
    # 課題 1-4
    return []

def translate(fastafile: str) -> List[str]:
    # 課題 1-5
    return []

if __name__ == "__main__":
    # filepath = "data/NT_113952.1.fasta"
    filepath = "data/ATGCCGT.fasta"
    # 課題 1-1
    print(base_count(filepath))
    # 課題 1-2
    print(gen_rev_comp_seq(filepath))
    # 課題 1-3
    print(calc_gc_content(filepath))
    # 課題 1-4
    print(search_motif(filepath, "ATG"))
    # 課題 1-5
    print(translate(filepath))
