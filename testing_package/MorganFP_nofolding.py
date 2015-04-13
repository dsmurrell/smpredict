# Oscar Mendez-Lucio 20130926
# Python 2.7.2
#
# Usage: phyton MorganFP_nofolding.py [MoleculesDataSet.sdf] [MoleculesExternalSet.sdf]
# In case a MoleculesExternalSet.sdf is not provided, the script generate the
# fingerprints for MoleculesDataSet.sdf
# Output:
# FPbits_Dictionary.txt - File containing SMILES for the fragment in each bit.
#                         It can be open with MarvinView (Chemaxon), openbabel
# BinaryFP.csv - File containing Morgan cicular fingerprints corresponding to the
#                molecular environment of each atom at certain radius (default=2).
#                The fragments and the number of bits depends on the molecules in
#                MoleculesDataSet.sdf.
#                
from __future__ import print_function
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit import DataStructs
from rdkit.Chem.Fingerprints import FingerprintMols
import sys

if len(sys.argv) == 2:
    sys.argv.append(sys.argv[1])

radius = 2 #It can be modified, radius=2 corresponds to ECFP4
DBcat = []
DataSet = sys.argv[1]
ExternalSet = sys.argv[2]


############################################################################
#
# Generation of Morgan circular fingerprints for MoleculesDataSet.sdf 
#
############################################################################

suppl = Chem.SDMolSupplier(DataSet)
for mol in suppl:
    for atom in range (0,mol.GetNumAtoms()):
        env = Chem.FindAtomEnvironmentOfRadiusN(mol,radius,atom)
        atomslist=set()
        for bidx in env:
            atomslist.add(mol.GetBondWithIdx(bidx).GetBeginAtomIdx())
            atomslist.add(mol.GetBondWithIdx(bidx).GetEndAtomIdx())
        smi = Chem.MolFragmentToSmiles(mol,atomsToUse=list(atomslist),bondsToUse=env,rootedAtAtom=atom)
        present = 0
        for i in range (0,len(DBcat)):
            if smi == DBcat[i]:
                present = 1
        if present == 0:
            DBcat.append(smi)
DBcat.sort()
print("Number of Bits:",len(DBcat))


############################################################################
#
# Generation of FPbits_Dictionary.txt
#
############################################################################

file = open("FPbits_Dictionary.txt", "w")
for i in range (0,len(DBcat)):
    print (DBcat[i],i+1,sep= "\t",file=file)
file.close()


############################################################################
#
# Generation of Morgan circular fingerprints for MoleculesExternalSet.sdf
# Generation of BinaryFP.csv
############################################################################

file= open("BinaryFP.csv", "w")
print ("\"Number of Bits:",len(DBcat),"\"",file=file)
suppl = Chem.SDMolSupplier(ExternalSet)
for mol in suppl:
    name = mol.GetProp("_Name")
    print ("\"",name,"\"", sep= "",end= ",",file=file)
    fp = ["0"] * len(DBcat)
    for atom in range (0,mol.GetNumAtoms()):
        env = Chem.FindAtomEnvironmentOfRadiusN(mol,radius,atom)
        atomslist=set()
        for bidx in env:
            atomslist.add(mol.GetBondWithIdx(bidx).GetBeginAtomIdx())
            atomslist.add(mol.GetBondWithIdx(bidx).GetEndAtomIdx())
        smi = Chem.MolFragmentToSmiles(mol,atomsToUse=list(atomslist),bondsToUse=env,rootedAtAtom=atom)
        for i in range (0,len(DBcat)):
            if smi == DBcat[i]:
                fp[i] = "1"
    print ("\"",','.join(fp),"\"",sep= "",file=file)
file.close()









