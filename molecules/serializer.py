from rest_framework import serializers
from .models import Molucule

class MoluculeSerializer(serializers.ModelSerializer):
    class Meta:
        model = Molucule
        fields = "__all__"