from django.shortcuts import render
from rdkit import Chem
from rdkit.Chem import AllChem
import py3Dmol
from django.http import JsonResponse
from django.views.decorators.csrf import csrf_exempt
import json            
import pubchempy as pcp

viewer_size = {'width': 1000, 'height': 600}  # Значення за замовчуванням

@csrf_exempt
def set_viewer_size(request):
    global viewer_size
    if request.method == 'POST':
        data = json.loads(request.body)
        viewer_size['width'] = int(data.get('width', 1000))
        viewer_size['height'] = int(data.get('height', 600))
    return JsonResponse({'status': 'success'})


def index(request):
    return render(request, 'index.html')

def about(request):
    return render(request, 'about.html')

def what_are_smiles(request):
    return render(request, 'what_are_smiles.html')


def molecule_view(request, name, smiles):
    global viewer_size
    width = viewer_size["width"]
    height = viewer_size["height"]
    
    print(width, height)
    
    mol = Chem.MolFromSmiles(smiles)
    mol = Chem.AddHs(mol)

    # Генеруємо 3D координати
    AllChem.EmbedMolecule(mol)
    block = Chem.MolToMolBlock(mol)

    # Візуалізація молекули
    viewer = py3Dmol.view(width=width, height=height)
    viewer.addModel(block, "mol")

    # Стиль молекули
    viewer.setStyle({'stick': {}})  # Відображаємо молекулу як стіки

    # Фоновий колір і зум
    viewer.setBackgroundColor('white')
    viewer.zoomTo()
    viewer.render()

    return render(request, 'molecule_view.html', {'mol_block': block, 'viewer': viewer, 'name': name})


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
