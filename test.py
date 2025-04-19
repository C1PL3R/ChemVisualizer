import pubchempy as pcp

def get_compounds_by_formula(formula):
    return pcp.get_compounds(formula, namespace='formula')

results = get_compounds_by_formula('RaSO4')

if results:
    compound = results[0]
    smiles = compound.isomeric_smiles
    molecular_formula = compound.molecular_formula
    molecular_weight = compound.molecular_weight
    exact_mass = compound.exact_mass
    iupac_name = compound.iupac_name
    isomeric_smiles = compound.isomeric_smiles
    aids = compound.aids if compound.aids else []
    sids = compound.sids if compound.sids else []
    synonyms = compound.synonyms if compound.synonyms else []

    print('Молекулярна формула:', molecular_formula)
    print('Молекулярна маса:', molecular_weight)
    print('SMILES:', isomeric_smiles)
    print("IUPAC-ім'я:", iupac_name)
    print('Точна маса:', exact_mass)
    print('AIDs:', aids[:3])  # показати перші 3 AID
    print('SIDs:', sids[:3])  # показати перші 3 SID
    print('Синоніми:', synonyms[:3])  # показати перші 3 назви
else:
    print("Молекулу не знайдено.")

