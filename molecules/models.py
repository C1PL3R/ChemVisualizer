from django.db import models


class Movie(models.Model):
    email = models.EmailField(max_length=20, null=True)
    password = models.CharField(max_length=30, null=True)
    email = models.EmailField(max_length=20, null=True)

    def __str__(self):
        return self.title