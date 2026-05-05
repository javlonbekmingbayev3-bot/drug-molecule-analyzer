from utils import molecular_weight, calculate_solubility, ro_5, display_formula, ADME
from rdkit import Chem

def main():
    print("===Welcome to Drug Molecule Analyzer program===")
    while True:
        molecule_structure = input("Enter a molecule in SMILES format: ")
        mol = Chem.MolFromSmiles(molecule_structure)
        if mol is not None:
            break
        else:
            print("Invalid SMILES string. Try again.")
            
    
    print("\nMolecular weight: "+ str(molecular_weight(mol)))
    print("Solubility: "+str(calculate_solubility(mol)))
    print("Lipinski Status: "+str(ro_5(mol)))
    print("===Optional Menu===")
    while True:
        print("1. Save 2D Visualization \n2. Check ADME (tPSA & Rotatable bonds) \n3. Exit")
        decision = int(input(": "))
        match decision:
            case 1:
                display_formula(mol)
            case 2:
                print(ADME(mol))
            case 3:
                break

if __name__ == "__main__":
    main()
