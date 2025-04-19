from django.urls import path
from django.conf import settings
from django.conf.urls.static import static
from django.conf.urls import handler400, handler403, handler404, handler500
from . import views
from rest_framework import routers
from .views import MoluculeAPIView


handler404 = 'your_project.views.custom_page_not_found'
handler500 = 'your_project.views.custom_error_view'

router = routers.DefaultRouter()
router.register(r'molecule-history', MoluculeAPIView, basename='molecule-history')

urlpatterns = [
    path('', views.index, name='index'),
    path('view_molecule/', views.molecule_view, name='view_molecule'),
    path('about/', views.about, name='about'),
    path('what_are_smiles/', views.what_are_smiles, name='what_are_smiles'),
    path('send-formula/', views.send_formula, name='send_formula'),
    path('send-name/', views.send_name, name='send-name')
] + static(settings.STATIC_URL, document_root=settings.STATIC_ROOT)