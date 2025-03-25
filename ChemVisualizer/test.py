from pubchempy import get_compounds

try:
    compounds = get_compounds('C8H8', 'formula')
    
    if compounds:
        for compound in compounds:
            print(compound.isomeric_smiles)
    else:
        print("Сполука не знайдена!")

except Exception as e:
    print(f"Помилка отримання даних з PubChem: {e}")