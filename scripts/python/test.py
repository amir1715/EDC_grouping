from glob import glob
import os
import pandas as pd
import numpy as np
import rdkit
from rdkit.Chem import AllChem
from rdkit.ML.Descriptors import MoleculeDescriptors
from rdkit.Chem import rdFingerprintGenerator
from rdkit.Chem import Descriptors
from rdkit import Chem

ligands = []
list_of_smiles = ['C1CCCC2C1CCCC2','C1CCCC2CCCCC12']

for i in list_of_smiles:
  mol = Chem.MolFromSmiles(i)
ligands.append(mol)

rdkit_df = pd.DataFrame(index=np.arange(len(ligands))) #empty dataframe to add the 200 rdkit features
fp_df = pd.DataFrame(index=np.arange(len(ligands)), columns=np.arange(124)) #dataframe to add the 124 bits of morgan fingerprints

for i in Descriptors._descList:
  rdkit_df[i[0]] = ' ' #create empty column with rdkit descriptor name#create empty column with rdkit descriptor name


morgan_fp_gen = Chem.rdFingerprintGenerator.GetMorganGenerator(includeChirality=True, radius=2, fpSize=124) #define morgan fingerprint generator
calc = MoleculeDescriptors.MolecularDescriptorCalculator([x[0] for x in Descriptors._descList])
#define the calculator for molecular descriptor

k = 0
for lig in ligands:
  desc = list(calc.CalcDescriptors(lig))
# print(desc)
fp = morgan_fp_gen.GetFingerprint(lig)
vector = np.array(fp)
print(len(vector))
rdkit_df.iloc[k, :] = desc #append the rdkit descriptors for each molecule in ligands
fp_df.iloc[k, :] = vector #append the bite vector for morgan fingerprint for each molecule in ligands
k += 1

