from glob import glob
import os
import pandas as pd
import numpy as np
import rdkit
from rdkit.Chem import AllChem
from rdkit.ML.Descriptors import MoleculeDescriptors
from rdkit.Chem import Descriptors
from rdkit import Chem

df=pd.read_csv('smiles.txt',sep='\t',header=None)
df.columns=['cas','smiles']

selected=df[~df.smiles.str.contains('not')]
list_of_smiles=selected['smiles'].tolist()
list_of_cas=selected['cas'].tolist()
ligands = []

for i in list_of_smiles:
  mol = Chem.MolFromSmiles(i)
  ligands.append(mol)

rdkit_df = pd.DataFrame(index=np.arange(len(ligands))) #empty dataframe to add the 200 rdkit features

for i in Descriptors._descList:
  rdkit_df[i[0]] = ' ' #create empty column with rdkit descriptor name#create empty column with rdkit descriptor name

calc = MoleculeDescriptors.MolecularDescriptorCalculator([x[0] for x in Descriptors._descList])
#define the calculator for molecular descriptor

k = 0
for lig in ligands:
  desc = list(calc.CalcDescriptors(lig))
  print(k)
  rdkit_df.iloc[k, :] = desc #append the rdkit descriptors for each molecule in ligands
  k += 1


rdkit_df['cas']=list_of_cas
rdkit_df.to_csv(r'RD_kit_res.csv',index=False,header=True)
