from django.shortcuts import render
from rdkit import Chem
from rdkit.Chem import AllChem, rdmolfiles
import py3Dmol
from django.http import JsonResponse, HttpResponse
from django.views.decorators.csrf import csrf_exempt
import json
import pubchempy as pcp
from .serializer import MoluculeSerializer
from .models import Molucule
from rest_framework import viewsets
from rest_framework.decorators import api_view
from rest_framework.response import Response


molecular_formula = None
molecular_weight = None
iupac_name = None
image_url = None
elements = None

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

    results = pcp.get_compounds(name, 'name')

    if results:
        compound = results[0]
        smiles = compound.isomeric_smiles
        molecular_formula = compound.molecular_formula
        molecular_weight = compound.molecular_weight
        exact_mass = compound.exact_mass
        iupac_name = compound.iupac_name
        aids = compound.aids if compound.aids else []
        sids = compound.sids if compound.sids else []        

    mol = Chem.MolFromSmiles(smiles)
    mol = Chem.AddHs(mol)

    # Генеруємо 3D координати
    AllChem.EmbedMolecule(mol)
    block = Chem.MolToMolBlock(mol)

    try:
        if name and name.strip():
            created = Molucule.objects.get_or_create(smiles=smiles, defaults={'name': str(name).capitalize()})
    except Exception as e:
        print(f"Error creating molecule: {e}")

    try:
        context = {
            'mol_block': block, 
            'name': name,
            'molecular_formula': molecular_formula,
            'molecular_weight': molecular_weight,
            'exact_mass': exact_mass,
            'iupac_name': iupac_name,
            'aids': aids[0],
            'sids': sids[0],
        }
    except IndexError:
        context = {
            'mol_block': block, 
            'name': name,
            'molecular_formula': molecular_formula,
            'molecular_weight': molecular_weight,
            'exact_mass': exact_mass,
            'iupac_name': iupac_name,
            'aids': None,
            'sids': None,
        }

    response = render(request, 'molecule_view.html', context)
    # response.delete_cookie('moleculeName', path='/')
    # response.delete_cookie('moleculeSmiles', path='/')
    return response


@api_view(['POST'])
def send_formula(request):
    formula = request.data.get('formula')
    results = pcp.get_compounds(formula, namespace='formula')
            
    if results:
        compound = results[0]
        smiles = compound.isomeric_smiles
        name = compound.synonyms[0]
        return JsonResponse({'status': 'success', 'name': name, 'smiles': smiles})
    else:
        return JsonResponse({'status': 'fail', 'error': 'Сполуку не знайдено'}, status=404)


@api_view(['POST'])
def send_name(request):
    name = request.data.get('name')
    results = pcp.get_compounds(name, 'name')

    if results:
        compound = results[0]
        smiles = compound.isomeric_smiles
        return JsonResponse({'status': 'success', 'name': name, 'smiles': smiles})
    else:
        return JsonResponse({'status': 'fail', 'error': 'Сполуку не знайдено'}, status=404)


def custom_page_not_found(request):
    return render(request, '404.html', status=404)


def custom_error_view(request):
    return render(request, '500.html', status=500)


class MoluculeAPIView(viewsets.ReadOnlyModelViewSet):
    queryset = Molucule.objects.all()
    serializer_class = MoluculeSerializer
