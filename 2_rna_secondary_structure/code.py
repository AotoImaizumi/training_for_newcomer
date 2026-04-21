from typing import List, Tuple, Union
import numpy.typing as npt
import numpy as np
from Bio import SeqIO
from Bio.PDB.PICIO import enumerate_atoms
from Bio.Seq import Seq

def supplement_RNA_base(base:str)->str:
    # RNAに限定したシンプルな相補塩基変換
    mapping = {"A": "U", "U": "A", "G": "C", "C": "G"}
    return mapping.get(base.upper(), "")

# fastaの配列は1つと仮定
def enumerate_pairs(fastafile: str) -> List[Tuple[int, int]]:
    # 課題 2-1
    with open(fastafile) as f:
        result = []
        for record in SeqIO.parse(f, "fasta"):
            # id = record.id
            # description = record.description
            seq = record.seq

            for i in range(len(seq)-1):
                supplement = supplement_RNA_base(seq[i])
                for j in range(i+1, len(seq)):
                    if(seq[j] == supplement):
                        result.append((i+1,j+1))
    return result

def enumerate_possible_pairs(fastafile: str, min_distance: int=4) -> List[Tuple[int, int]]:
    # 課題 2-2
    # 間に4塩基以上挟まっているという認識でいいのか？
    return [tp for tp in enumerate_pairs(fastafile) if tp[1]-tp[0]>4]


def enumerate_continuous_pairs(fastafile: str, min_distance: int=4, min_length: int=2) -> List[Tuple[int, int, int]]:
    # 課題 2-3
    enumerate_possible_pairs(fastafile, min_distance)
    return []

def create_dotbracket_notation(fastafile: str, min_distance: int=4, min_length: int=2) -> str:
    # 課題 2-4
    return ""

if __name__ == "__main__":
    filepath = "data/AUCGCCAU.fasta"
    # 課題 2-1
    print(enumerate_pairs(filepath))
    # 課題 2-2
    print(enumerate_possible_pairs(filepath))
    # 課題 2-3
    print(enumerate_continuous_pairs(filepath, 2))
    # 課題 2-4
    print(create_dotbracket_notation(filepath, 2))


