from utils import molecular_weight, calculate_solubility, ro_5, display_formula, ADME

def main():
    print("===Welcome to Drug Molecule Analyzer program===")
    molecule_structure = input("Enter a molecule in SMILES format: ")
    
    print("Molecule weight: "+ str(molecular_weight(molecule_structure)))
    print("Solubility: "+str(calculate_solubility(molecule_structure)))
    print("Lipinski Status: "+str(ro_5(molecule_structure)))
    print("===Optional Menu===")
    while True:
        print("1. Save 2D Visualization \n2. Check ADME (tPSA & Rotatable bonds) \n3. Exit")
        decision = int(input(": "))
        match decision:
            case 1:
                print(display_formula(molecule_structure))
            case 2:
                print(ADME(molecule_structure))
            case 3:
                break

if __name__ == "__main__":
    main()
