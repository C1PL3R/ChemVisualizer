from django.urls import path
from . import views

urlpatterns = [
    path('', views.index, name='index'),
    path('molecule-smiles/', views.get_molecule_smiles, name='get_molecule_smiles'),
    path('view_molecule/<path:name>/<path:smiles>/', views.molecule_view, name='view_molecule'),
    path('set_viewer_size/', views.set_viewer_size, name='set_viewer_size'),
    path('about/', views.about, name='about'),
    path('what_are_smiles/', views.what_are_smiles, name='what_are_smiles')
]
