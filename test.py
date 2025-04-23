from rdkit import Chem
from rdkit.Chem import AllChem
from stl import mesh
import numpy as np
from math import pi, sin, cos

# === Налаштування ===
SDF_FILE = "molecule.sdf"
STL_FILE = "molecule.stl"
ATOM_RADIUS = 0.3
BOND_RADIUS = 0.07
FACES = 20

def create_sphere(center, radius, faces=FACES):
    vertices = []
    triangles = []
    for i in range(faces):
        lat0 = pi * (-0.5 + float(i) / faces)
        z0 = radius * sin(lat0)
        zr0 = radius * cos(lat0)

        lat1 = pi * (-0.5 + float(i + 1) / faces)
        z1 = radius * sin(lat1)
        zr1 = radius * cos(lat1)

        for j in range(faces):
            lng0 = 2 * pi * float(j) / faces
            lng1 = 2 * pi * float(j + 1) / faces

            x0 = cos(lng0)
            y0 = sin(lng0)
            x1 = cos(lng1)
            y1 = sin(lng1)

            p1 = [center[0] + zr0 * x0, center[1] + zr0 * y0, center[2] + z0]
            p2 = [center[0] + zr1 * x0, center[1] + zr1 * y0, center[2] + z1]
            p3 = [center[0] + zr1 * x1, center[1] + zr1 * y1, center[2] + z1]
            p4 = [center[0] + zr0 * x1, center[1] + zr0 * y1, center[2] + z0]

            vertices.extend([p1, p2, p3, p4])
            base = len(vertices) - 4
            triangles.append([base, base + 1, base + 2])
            triangles.append([base, base + 2, base + 3])
    return vertices, triangles

# === Створення циліндра ===
def create_cylinder(start, end, radius=BOND_RADIUS, segments=12):
    start = np.array(start)
    end = np.array(end)
    axis = end - start
    length = np.linalg.norm(axis)
    if length == 0:
        return [], []
    axis /= length
    not_axis = np.array([1, 0, 0]) if abs(axis[0]) < 0.99 else np.array([0, 1, 0])
    v = np.cross(axis, not_axis)
    v /= np.linalg.norm(v)
    w = np.cross(axis, v)

    angle_step = 2 * pi / segments
    vertices = []
    triangles = []

    for i in range(segments):
        angle = i * angle_step
        next_angle = (i + 1) * angle_step

        p1 = start + radius * (cos(angle) * v + sin(angle) * w)
        p2 = start + radius * (cos(next_angle) * v + sin(next_angle) * w)
        p3 = p1 + axis * length
        p4 = p2 + axis * length

        base = len(vertices)
        vertices.extend([p1, p2, p3, p4])
        triangles.append([base, base + 1, base + 2])
        triangles.append([base + 1, base + 3, base + 2])
    return vertices, triangles

# === Завантаження молекули ===
mol = Chem.MolFromMolFile(SDF_FILE, removeHs=False)
if not mol:
    raise ValueError("Не вдалося завантажити молекулу з SDF.")

if not mol.GetConformer().Is3D():
    AllChem.EmbedMolecule(mol)
    AllChem.UFFOptimizeMolecule(mol)

conf = mol.GetConformer()
all_vertices = []
all_triangles = []

# === Сфери для атомів ===
for atom in mol.GetAtoms():
    pos = conf.GetAtomPosition(atom.GetIdx())
    center = [pos.x, pos.y, pos.z]
    vertices, triangles = create_sphere(center, ATOM_RADIUS)
    offset = len(all_vertices)
    all_vertices.extend(vertices)
    all_triangles.extend([[a + offset for a in tri] for tri in triangles])

# === Циліндри для зв’язків ===
for bond in mol.GetBonds():
    start = conf.GetAtomPosition(bond.GetBeginAtomIdx())
    end = conf.GetAtomPosition(bond.GetEndAtomIdx())
    start_coords = [start.x, start.y, start.z]
    end_coords = [end.x, end.y, end.z]
    vertices, triangles = create_cylinder(start_coords, end_coords)
    offset = len(all_vertices)
    all_vertices.extend(vertices)
    all_triangles.extend([[a + offset for a in tri] for tri in triangles])

# === Експорт STL ===
all_vertices = np.array(all_vertices)
all_triangles = np.array(all_triangles)

mol_mesh = mesh.Mesh(np.zeros(all_triangles.shape[0], dtype=mesh.Mesh.dtype))
for i, f in enumerate(all_triangles):
    for j in range(3):
        mol_mesh.vectors[i][j] = all_vertices[f[j]]

mol_mesh.save(STL_FILE)
print(f"✅ STL файл збережено: {STL_FILE}")
