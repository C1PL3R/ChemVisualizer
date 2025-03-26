from django.shortcuts import render
from rdkit import Chem
from rdkit.Chem import AllChem
import py3Dmol
from django.http import JsonResponse
from django.views.decorators.csrf import csrf_exempt
import json            
import pubchempy as pcp            
import urllib.parse


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
            decoded_str = urllib.parse.unquote(name)
            
            def get_compounds_by_formula(formula):
                compounds = pcp.get_compounds(formula, 'name')
                return compounds

            results = get_compounds_by_formula(decoded_str)
            
            smiles = None  # Ініціалізація змінної перед циклом
            
            if results:
                for compound in results:
                    if compound.isomeric_smiles:
                        smiles = compound.isomeric_smiles
                        break  # Беремо перший знайдений SMILES і виходимо з циклу
            if smiles:
                return JsonResponse({'status': 'success', 'smiles': smiles, 'name': decoded_str})
            else:
                return JsonResponse({'status': 'fail', 'error': 'SMILES not found for given name'}, status=404)
        
        except json.JSONDecodeError:
            return JsonResponse({'status': 'fail', 'error': 'Invalid JSON'}, status=400)
        except Exception as e:
            print(e)
            return JsonResponse({'status': 'fail', 'error': str(e)}, status=500)

    return JsonResponse({'status': 'fail', 'error': 'Invalid HTTP method'}, status=405)


def custom_page_not_found(request, exception):
    return render(request, '404.html', status=404)

def custom_error_view(request):
    return render(request, '500.html', status=500)