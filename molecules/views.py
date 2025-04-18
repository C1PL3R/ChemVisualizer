from django.shortcuts import render
from rdkit import Chem
from rdkit.Chem import AllChem, rdmolfiles
import py3Dmol
from django.http import JsonResponse
from django.views.decorators.csrf import csrf_exempt
import json
import pubchempy as pcp
from .serializer import MoluculeSerializer
from .models import Molucule
from rest_framework import viewsets
from django.contrib import messages


def index(request):
    molucules = Molucule.objects.all()
    return render(request, 'index.html', {'molecules_history': molucules})


def about(request):
    return render(request, 'about.html')


def what_are_smiles(request):
    return render(request, 'what_are_smiles.html')


def molecule_view(request):
    name = request.COOKIES.get('moleculeName')
    smiles = request.COOKIES.get('moleculeSmiles')


    if smiles is None or smiles == "":
        def get_compounds_by_formula(formula):
            compounds = pcp.get_compounds(formula, 'name')
            return compounds

        results = get_compounds_by_formula(name)

        for compound in results:
            smiles = compound.isomeric_smiles

    mol = Chem.MolFromSmiles(smiles)
    mol = Chem.AddHs(mol)

    # Генеруємо 3D координати
    AllChem.EmbedMolecule(mol)
    block = Chem.MolToMolBlock(mol)
    
    try:
        if name and name.strip():
            Molucule.objects.get_or_create(name=str(name).capitalize, smiles=smiles)
    except Exception as e:
        print(f"Error creating molecule: {e}")

    response = render(request, 'molecule_view.html', {
                      'mol_block': block, 'name': name})
    response.delete_cookie('moleculeName', path='/')
    response.delete_cookie('moleculeSmiles', path='/')
    return response


@csrf_exempt
def check_smiles_code(request):
    if request.method == 'POST':
        try:
            data = json.loads(request.body)
            smiles = data.get("smiles")

            if not smiles or smiles.strip() == "":
                return JsonResponse({'status': 'fail', 'error': 'SMILES code is empty!'}, status=400)

            def is_valid_smiles(smiles):
                mol = Chem.MolFromSmiles(smiles)
                return mol is not None

            if is_valid_smiles(smiles):
                return JsonResponse({'status': 'success'})
            else:
                return JsonResponse({'status': 'fail', 'error': 'SMILES code is enter incorrectly!'}, status=404)

        except json.JSONDecodeError:
            return JsonResponse({'status': 'fail', 'error': 'Invalid JSON'}, status=400)
        except Exception as e:
            return JsonResponse({'status': 'fail', 'error': str(e)}, status=500)

    return JsonResponse({'status': 'fail', 'error': 'Invalid HTTP method'}, status=405)


def custom_page_not_found(request):
    return render(request, '404.html', status=404)


def custom_error_view(request):
    return render(request, '500.html', status=500)


class MoluculeAPIView(viewsets.ModelViewSet):
    queryset = Molucule.objects.all()
    serializer_class = MoluculeSerializer
