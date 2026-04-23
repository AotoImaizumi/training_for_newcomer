from typing import List, Union
import numpy.typing as npt
import numpy as np
import pandas as pd
from rdkit import Chem
from rdkit.Chem import Draw
from sklearn.model_selection import train_test_split
from sklearn.ensemble import RandomForestRegressor

#4-0
"""
細胞膜が特定の物質をどれだけ透過させやすいかという指標
創薬の観点での役割；物質を細胞内に効率的に浸透させたい
experimentalに調べるには時間とコストがかかりすぎるので、スクリーニングしたい
中分子は大きさ的に膜を通りにくいため最適な構造を事前に調べて効率化したい
"""

def draw_molecule(csvfile: str, molecule:str) -> None:
    # 課題 4-1
    # RDKitで描画可
    df = pd.read_csv(csvfile)
    smiles = df.query("`Compound ID` == @molecule")['SMILES'].item()
    mol = Chem.MolFromSmiles(smiles)
    Draw.MolToFile(mol, filename=f"{molecule}.svg")
    return

def create_2d_descriptors(smiles: str) -> Union[npt.NDArray[np.float64], List[float]]:
    # 課題 4-2
    return []

def predict_logpapp(csvfile: str) -> Union[npt.NDArray[np.float64], pd.Series, List[float]]:
    # 課題 4-3
    np.random.seed(0) # 出力を固定するためにseedを指定
    rfr = RandomForestRegressor(random_state=0) # 出力を固定するためにrandom_stateを指定

    return 0.0

def grid_search(csvfile: str) -> float:
    # 課題 4-4
    # こちらも出力を固定するためにseedやrandom_stateを指定すること
    np.random.seed(0)

    # # Xは説明変数、yは目的変数
    # X_train, X_test, y_train, y_test = train_test_split(X, y, train_size=700, random_state=0)
    # 

    return 0.0

if __name__ == "__main__":
    smiles = "C(=O)(c1ccc(OCCCCCC)cc1)CCNc1cc(Cl)ccc1"
    filepath = "data/fukunishi_data.csv"
    # 課題 4-1
    draw_molecule(filepath, "CHEMBL540227")
    # 課題 4-2
    print(create_2d_descriptors(smiles))
    # 課題 4-3
    print(predict_logpapp(filepath))
    # 課題 4-4
    print(grid_search(filepath))
