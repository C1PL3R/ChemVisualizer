from django.shortcuts import render, redirect
from rdkit import Chem
from rdkit.Chem import AllChem

from django.http import JsonResponse, HttpResponse

import pubchempy as pcp
from pubchempy import BadRequestError
from rest_framework import viewsets
from rest_framework.decorators import api_view

from .serializer import MoluculeSerializer
from .models import Molucule
from .converter import convert_sdf_to_stl_bytes

import os


molecular_formula = None
molecular_weight = None
iupac_name = None
image_url = None
elements = None

def converter_to_sdf(request):
    if request.method == "POST" and request.FILES.get("sdf_file"):
        sdf_file = request.FILES["sdf_file"]
        sdf_data = sdf_file.read().decode("utf-8")
        
        try:
            stl_data = convert_sdf_to_stl_bytes(sdf_data)

            # Отримуємо ім'я файлу без розширення
            filename = os.path.splitext(sdf_file.name)[0]

            return HttpResponse(stl_data, content_type="application/sla", headers={
                'Content-Disposition': f'attachment; filename="{filename}.stl"'
            })
        except ValueError as e:
            return HttpResponse(f"❌ Помилка: {str(e)}", status=400)

    return render(request, 'converter.html')


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
        synonyms = compound.synonyms if compound.synonyms else []     

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
            'synonyms': f"{synonyms[0]}, {synonyms[1]}, {synonyms[3]}",
        }
    except IndexError:
        context = {
            'mol_block': block, 
            'name': name,
            'molecular_formula': molecular_formula,
            'molecular_weight': molecular_weight,
            'exact_mass': exact_mass,
            'iupac_name': iupac_name,
            'synonyms': f"{synonyms[0]}, {synonyms[1]}, {synonyms[3]}",
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
    try:
        if results:
            compound = results[0]
            smiles = compound.isomeric_smiles
            name = compound.synonyms[0]
            return JsonResponse({'status': 'success', 'name': name, 'smiles': smiles}, status=200)
        else:
            return JsonResponse({'status': 'fail', 'error': 'Сполуку не знайдено'})    
    except BadRequestError:
        return JsonResponse({'status': 'fail', 'error': 'Некоректна формула'})


@api_view(['POST'])
def send_name(request):
    name = request.data.get('name')
    results = pcp.get_compounds(name, 'name')

    if results:
        compound = results[0]
        smiles = compound.isomeric_smiles
        return JsonResponse({'status': 'success', 'name': name, 'smiles': smiles}, status=200)
    else:
        return JsonResponse({'status': 'fail', 'error': 'Сполуку не знайдено'})


def custom_page_not_found(request):
    return render(request, '404.html', status=404)


def custom_error_view(request):
    return render(request, '500.html', status=500)


def visualize_sdf(request):
    if request.method == 'POST' and request.FILES.get('sdf_file'):
        sdf_file = request.FILES['sdf_file']
        sdf_data = sdf_file.read().decode('utf-8')
        request.session['sdf_data'] = sdf_data  # тимчасово в сесію
        return redirect('visualize_sdf')  # після POST → GET

    sdf_data = request.session.pop('sdf_data', None)
    return render(request, 'visualize_sdf.html', {'sdf_data': sdf_data})


class MoluculeAPIView(viewsets.ReadOnlyModelViewSet):
    queryset = Molucule.objects.all()
    serializer_class = MoluculeSerializer
