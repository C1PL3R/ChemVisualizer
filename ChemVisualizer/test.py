import pubchempy as pcp

def get_compounds_by_formula(formula):
    compounds = pcp.get_compounds(formula, 'name')
    return compounds

# Отримуємо список з'єднань для формули C6H12O6 (може бути кілька ізомерів)
results = get_compounds_by_formula("Glucose")

# Виводимо інформацію
if results:
    for compound in results:
        print(f"CID: {compound.cid}")
        print(f"Назва: {compound.iupac_name}")
        print(f"Формула: {compound.molecular_formula}")
        print(f"SMILES: {compound.isomeric_smiles}")
        print(f"InChI: {compound.inchi}\n")
else:
    print("Нічого не знайдено в PubChem.")
