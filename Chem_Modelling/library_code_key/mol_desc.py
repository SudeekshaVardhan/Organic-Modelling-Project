# import libraries necessary
import pandas as pd
import warnings

from rdkit import Chem
from rdkit.Chem import Draw
from rdkit.Chem import Descriptors
from rdkit.Chem import AllChem, PandasTools
from rdkit.Chem import Descriptors

warnings.filterwarnings("ignore")

# 
file_path = r'C:\Users\Sudee\Documents\coding_projects\Chem_Modelling\library_code_key\Heteroaromatics.xlsx'
df = pd.read_excel(file_path)

PandasTools.AddMoleculeColumnToFrame(df, 'Smiles', 'mol')
df

mol_list = []

for smile in df['Smiles']:
    mol = Chem.MolFromSmiles(smile)
    mol = Chem.AddHs(mol)
    mol_list = list.append(mol)
# creating molecular object from SMILES
