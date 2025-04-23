from django.db import models

    
class Molucule(models.Model):
    name = models.CharField(max_length=50, null=False, unique=True)
    smiles = models.TextField(null=False, unique=True)
    created_at = models.DateTimeField(auto_now_add=True)
    
    def __str__(self):
        return self.name