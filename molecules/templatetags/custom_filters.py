from django import template
from django.utils.safestring import mark_safe
import re

register = template.Library()

@register.filter
def format_formula(formula):
    """
    Перетворює C6H12O6 → C<sub>6</sub>H<sub>12</sub>O<sub>6</sub>
    """
    if not isinstance(formula, str):
        return formula
    formatted = re.sub(r'(\d+)', r'<sub>\1</sub>', formula)
    return mark_safe(formatted)