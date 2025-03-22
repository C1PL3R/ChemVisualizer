from django.shortcuts import render
from rdkit import Chem
from rdkit.Chem import AllChem
import py3Dmol
from django.http import JsonResponse
from django.views.decorators.csrf import csrf_exempt
import json            
import pubchempy as pcp

# Create your views here.
def index(request):
    return render(request, 'index.html')


def molecule_view(request, smiles):
    # Конвертуємо SMILES в молекулу
    mol = Chem.MolFromSmiles(smiles)
    mol = Chem.AddHs(mol)

    # Генеруємо 3D координати
    AllChem.EmbedMolecule(mol)
    block = Chem.MolToMolBlock(mol)

    # Візуалізація молекули
    viewer = py3Dmol.view(width=1000, height=600)
    viewer.addModel(block, "mol")

    # Стиль молекули
    viewer.setStyle({'stick': {}})  # Відображаємо молекулу як стіки

    # Додаємо підписи атомів
    for atom in mol.GetAtoms():
        pos = mol.GetConformer().GetAtomPosition(atom.GetIdx())
        viewer.addLabel(atom.GetSymbol(), {
            "position": {'x': pos.x, 'y': pos.y, 'z': pos.z},
            "fontSize": 12,
            "fontColor": "black"
        })

    # Фоновий колір і зум
    viewer.setBackgroundColor('white')
    viewer.zoomTo()
    viewer.render()

    return render(request, 'molecule_view.html', {'mol_block': block, 'viewer': viewer})


@csrf_exempt
def get_molecule_smiles(request):
    if request.method == 'POST':
        try:
            data = json.loads(request.body)
            name = data.get("name")

            def get_compounds_by_formula(formula):
                compounds = pcp.get_compounds(formula, 'name')
                return compounds

            results = get_compounds_by_formula(name)
            
            if results:
                for compound in results:
                    smiles = compound.isomeric_smiles
            
            return JsonResponse({'status': 'success', 'smiles': smiles, 'name': name})
        except json.JSONDecodeError:
            return JsonResponse({'status': 'fail', 'error': 'Invalid JSON'}, status=400)
        except Exception as e:
            return JsonResponse({'status': 'fail', 'error': str(e)}, status=500)
    return JsonResponse({'status': 'fail', 'error': 'Invalid HTTP method'}, status=405)
