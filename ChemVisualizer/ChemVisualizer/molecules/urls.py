from django.urls import path
from . import views

urlpatterns = [
    path('', views.index, name='index'),
    path('molecule-smiles/', views.get_molecule_smiles, name='get_molecule_smiles'),
    path('view_molecule/<str:smiles>/', views.molecule_view, name='view_molecule')
]
