from rdkit import Chem
from rdkit import DataStructs
from rdkit.Chem import AllChem

# Define SMILES strings
smiles1 = "CC(C)C=CCCCCC(=O)NCc1ccc(c(c1)OC)O"
smiles2 = "COC1=C(C=CC(=C1)C=O)O"

# Load molecules from SMILES strings
mol1 = Chem.MolFromSmiles(smiles1)
mol2 = Chem.MolFromSmiles(smiles2)

# Aromatize molecules (if needed)
mol1 = Chem.Mol(mol1)
mol2 = Chem.Mol(mol2)
Chem.SanitizeMol(mol1)
Chem.SanitizeMol(mol2)
Chem.Kekulize(mol1)
Chem.Kekulize(mol2)

# Calculate similarity between "Morgan" fingerprints
print("Morgan fingerprints:")
fp1 = AllChem.GetMorganFingerprintAsBitVect(mol1, 2)  # 2 is radius
fp2 = AllChem.GetMorganFingerprintAsBitVect(mol2, 2)

print("  Tanimoto: %s" % (DataStructs.TanimotoSimilarity(fp1, fp2)))

# Calculate similarity between "Substructure" fingerprints
print("Substructure fingerprints:")
sub_fp1 = Chem.PatternFingerprint(mol1)
sub_fp2 = Chem.PatternFingerprint(mol2)

print("  Tanimoto: %s" % (DataStructs.TanimotoSimilarity(sub_fp1, sub_fp2)))


# m1 = indigo.loadMolecule("CC(C)C=CCCCCC(=O)NCc1ccc(c(c1)OC)O")
# m2 = indigo.loadMolecule("COC1=C(C=CC(=C1)C=O)O")
# # Aromatize molecules because second molecule is not in aromatic form
# m1.aromatize()
# m2.aromatize()
#
# # Calculate similarity between "similarity" fingerprints
# print("Similarity fingerprints:");
# fp1 = m1.fingerprint("sim");
# fp2 = m2.fingerprint("sim");
#
# print("  Tanimoto: %s" % (indigo.similarity(fp1, fp2, "tanimoto")));
# print("  Tversky: %s" % (indigo.similarity(fp1, fp2, "tversky")));
#
# # Calculate similarity between "substructure" fingerprints
# print("Substructure fingerprints:");
# fp1 = m1.fingerprint("sub");
# fp2 = m2.fingerprint("sub");
#
# print("  Tanimoto: %s" % (indigo.similarity(fp1, fp2, "tanimoto")));
# print("  Tversky: %s" % (indigo.similarity(fp1, fp2, "tversky")));