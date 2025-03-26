from django.urls import path
from . import views
from django.conf import settings
from django.conf.urls.static import static
from django.conf.urls import handler400, handler403, handler404, handler500

handler404 = 'your_project.views.custom_page_not_found'
handler500 = 'your_project.views.custom_error_view'

urlpatterns = [
    path('', views.index, name='index'),
    path('molecule-smiles/', views.get_molecule_smiles, name='get_molecule_smiles'),
    path('view_molecule/<path:name>/<path:smiles>/', views.molecule_view, name='view_molecule'),
    path('set_viewer_size/', views.set_viewer_size, name='set_viewer_size'),
    path('about/', views.about, name='about'),
    path('what_are_smiles/', views.what_are_smiles, name='what_are_smiles'),
] + static(settings.STATIC_URL, document_root=settings.STATIC_ROOT)