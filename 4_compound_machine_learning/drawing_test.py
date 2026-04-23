from rdkit import Chem
from rdkit.Chem import Draw


mol = Chem.MolFromSmiles('C(=O)(c1ccc(OCCCCCC)cc1)CCNc1cc(Cl)ccc1	')
Draw.MolToFile( mol , filename='test.png')
