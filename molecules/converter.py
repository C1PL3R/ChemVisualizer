import tempfile
from django.http import HttpResponse
from rdkit import Chem
from rdkit.Chem import AllChem
from stl import mesh
import numpy as np
from math import pi, sin, cos

def convert_sdf_to_stl_bytes(sdf_data: str) -> bytes:
    mol = Chem.MolFromMolBlock(sdf_data, removeHs=False)
    if mol is None:
        raise ValueError("Не вдалося зчитати файл SDF")

    if not mol.GetConformer().Is3D():
        AllChem.EmbedMolecule(mol)
        AllChem.UFFOptimizeMolecule(mol)

    conf = mol.GetConformer()
    all_vertices = []
    all_triangles = []

    def create_sphere(center, radius, faces=20):
        vertices, triangles = [], []
        for i in range(faces):
            lat0 = pi * (-0.5 + i / faces)
            z0 = radius * sin(lat0)
            zr0 = radius * cos(lat0)
            lat1 = pi * (-0.5 + (i + 1) / faces)
            z1 = radius * sin(lat1)
            zr1 = radius * cos(lat1)
            for j in range(faces):
                lng0 = 2 * pi * j / faces
                lng1 = 2 * pi * (j + 1) / faces
                x0, y0 = cos(lng0), sin(lng0)
                x1, y1 = cos(lng1), sin(lng1)
                p1 = [center[0] + zr0 * x0, center[1] + zr0 * y0, center[2] + z0]
                p2 = [center[0] + zr1 * x0, center[1] + zr1 * y0, center[2] + z1]
                p3 = [center[0] + zr1 * x1, center[1] + zr1 * y1, center[2] + z1]
                p4 = [center[0] + zr0 * x1, center[1] + zr0 * y1, center[2] + z0]
                vertices.extend([p1, p2, p3, p4])
                base = len(vertices) - 4
                triangles.append([base, base + 1, base + 2])
                triangles.append([base, base + 2, base + 3])
        return vertices, triangles

    def create_cylinder(start, end, radius=0.1, segments=12):
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
        vertices, triangles = [], []
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

    # Створення сфери для атомів
    for atom in mol.GetAtoms():
        pos = conf.GetAtomPosition(atom.GetIdx())
        center = [pos.x, pos.y, pos.z]
        v, t = create_sphere(center, radius=0.3)
        offset = len(all_vertices)
        all_vertices.extend(v)
        all_triangles.extend([[i + offset for i in tri] for tri in t])

    # Створення циліндрів для зв’язків
    for bond in mol.GetBonds():
        start = conf.GetAtomPosition(bond.GetBeginAtomIdx())
        end = conf.GetAtomPosition(bond.GetEndAtomIdx())
        v, t = create_cylinder([start.x, start.y, start.z], [end.x, end.y, end.z])
        offset = len(all_vertices)
        all_vertices.extend(v)
        all_triangles.extend([[i + offset for i in tri] for tri in t])

    all_vertices = np.array(all_vertices)
    all_triangles = np.array(all_triangles)

    # Створення STL в тимчасовий файл
    with tempfile.NamedTemporaryFile(suffix=".stl", delete=False) as tmp:
        mol_mesh = mesh.Mesh(np.zeros(all_triangles.shape[0], dtype=mesh.Mesh.dtype))
        for i, f in enumerate(all_triangles):
            for j in range(3):
                mol_mesh.vectors[i][j] = all_vertices[f[j]]
        mol_mesh.save(tmp.name)  # Зберігаємо STL в тимчасовий файл
        tmp.seek(0)  # Переміщаємо в початок файлу
        stl_bytes = tmp.read()  # Читаємо в пам'ять як bytes

    return stl_bytes