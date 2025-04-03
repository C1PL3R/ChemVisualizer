import sys
import trimesh
import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem.rdchem import Mol

def mol_to_mesh(mol: Mol):
    mol = Chem.AddHs(mol)
    AllChem.EmbedMolecule(mol, AllChem.ETKDG())
    conf = mol.GetConformer()
    
    vertices = [list(conf.GetAtomPosition(i)) for i in range(mol.GetNumAtoms())]
    faces = []
    
    for bond in mol.GetBonds():
        i, j = bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()
        faces.append([i, j, (i + j) % mol.GetNumAtoms()])
    
    return trimesh.Trimesh(vertices=np.array(vertices), faces=np.array(faces))

def convert_mol_to_stl(mol_file: str, stl_file: str):
    mol = Chem.MolFromMolFile(mol_file)
    if not mol:
        raise ValueError("Не вдалося завантажити .mol файл")
    
    mesh = mol_to_mesh(mol)
    mesh.export(stl_file)
    print(f"Файл {stl_file} успішно створено!")

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Використання: python script.py input.mol output.stl")
        sys.exit(1)
    
    input_mol = sys.argv[1]
    output_stl = sys.argv[2]
    
    convert_mol_to_stl(input_mol, output_stl)
