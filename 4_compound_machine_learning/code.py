from cmath import sqrt
from typing import List, Union
import numpy.typing as npt
import numpy as np
import pandas as pd
from rdkit import Chem
from rdkit.Chem import Draw
from rdkit.Chem import Descriptors
from rdkit.ML.Descriptors import MoleculeDescriptors
from rdkit.ML.Descriptors.MoleculeDescriptors import MolecularDescriptorCalculator
from sklearn.metrics import mean_squared_error
from sklearn.model_selection import train_test_split
from sklearn.ensemble import RandomForestRegressor
from sklearn.model_selection import cross_val_score
from sklearn.model_selection import GridSearchCV

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
    # 1つのSMILESに対して2d記述子を作り出す
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None

    descriptor_names = [d[0] for d in Descriptors._descList]

    calculator = MoleculeDescriptors.MolecularDescriptorCalculator(descriptor_names)

    descriptor_values = calculator.CalcDescriptors(mol)
    return np.array(descriptor_values)

def predict_logpapp(csvfile: str) -> Union[npt.NDArray[np.float64], pd.Series, List[float]]:
    # 課題 4-3
    np.random.seed(0) # 出力を固定するためにseedを指定
    rfr = RandomForestRegressor(random_state=0) # 出力を固定するためにrandom_stateを指定
    pass

    # Xは説明変数、yは目的変数
    df = pd.read_csv(csvfile)
    X = np.array([create_2d_descriptors(sml) for sml in df["SMILES"]])
    y = df["LogP app"]
    X_train, X_test, y_train, y_test, id_train, id_test = train_test_split(X, y, df["No."], train_size=700, random_state=0)

    rfr.fit(X_train, y_train)
    y_pred = rfr.predict(X_test)

    result_df = pd.DataFrame({
        "No.": id_test,
        "Predicted_LogP": y_pred
    })
    return result_df.to_string(index=False)


def grid_search(csvfile: str) -> float:
    # 課題 4-4
    # こちらも出力を固定するためにseedやrandom_stateを指定すること
    np.random.seed(0)


    # # Xは説明変数、yは目的変数
    # X_train, X_test, y_train, y_test = train_test_split(X, y, train_size=700, random_state=0)
    df = pd.read_csv(csvfile)
    X = np.array([create_2d_descriptors(sml) for sml in df["SMILES"]])
    y = df["LogP app"]
    X_train, X_test, y_train, y_test, id_train, id_test = train_test_split(X, y, df["No."], train_size=700, random_state=0)

    param_grid = {
        'n_estimators': [100, 200, 400],
        'max_depth': [5, 10, 15]
    }

    # GridSearchCV のインスタンスを作成
    # cv=4: 4-fold交差検証, n_jobs=-1: 全コア使用
    grid_search = GridSearchCV(
        RandomForestRegressor(random_state=0),
        param_grid,
        cv=4,
        scoring='neg_mean_squared_error',
        n_jobs=-1
    )
    grid_search.fit(X_train, y_train)

    # 結果の確認
    # cvに用いた525件で学習し、175件のtest dataでのRMSEが出る
    print(f"最良の組み合わせ: {grid_search.best_params_}")
    print(f"その時のRMSE: {np.sqrt(-grid_search.best_score_):.3f}")

    # そのまま予測に使える（内部で最良モデルが再学習済み）
    # 最良モデルで700件のデータで学習し、94件のtest dataでのRMSEが出る→訓練件数が多いのでRMSEが小さく出る
    y_pred = grid_search.predict(X_test)
    rmse = sqrt(mean_squared_error(y_test, y_pred))
    return rmse

if __name__ == "__main__":
    np.set_printoptions(suppress=True, precision=17)
    smiles = "C(=O)(c1ccc(OCCCCCC)cc1)CCNc1cc(Cl)ccc1"
    filepath = "data/fukunishi_data.csv"
    # 課題 4-1
    #draw_molecule(filepath, "CHEMBL540227")
    # 課題 4-2
    #print("4-2")
    #print(create_2d_descriptors(smiles))
    # 課題 4-3
    #print("4-3")
    #print(predict_logpapp(filepath))
    # 課題 4-4
    print("4-4")
    print(grid_search(filepath))
