from rdkit.Chem import Descriptors, Lipinski, Draw
from rdkit import Chem

def molecular_weight(molecule):
    """Calculates the molecular weight of the given molecule in SMILES format."""
    mol = Chem.MolFromSmiles(molecule)
    m_w = Descriptors.MolWt(mol)
    return m_w

def calculate_solubility(molecule):
    """Calculates the solubility of the molecule."""
    mol = Chem.MolFromSmiles(molecule)
    logP = Descriptors.MolLogP(mol)
    if logP<1:
        return 'very soluble'
    elif logP> 1 and logP<3:
        return 'moderate soluble'
    else:
        return 'poorly soluble'
    
def ro_5(molecule):
    """Applies Lipinski's rule of five to test whether the molecule is orally bioavailable to the human body."""
    mol = Chem.MolFromSmiles(molecule)
    hbd_count = Lipinski.NumHDonors(mol)
    hba_count = Lipinski.NumHAcceptors(mol)
    m_w = molecular_weight(molecule)
    LogP = Descriptors.MolLogP(mol)

    violations = 0

    if m_w > 500:
        violations+=1
    if hba_count > 10:
        violations+=1
    if hbd_count > 5:
        violations+=1
    if LogP > 5:
        violations+=1

    return True if violations <= 1 else False
    
def display_formula(molecule):
    """Creates an image that shows displayed formula of the molecule."""
    mol = Chem.MolFromSmiles(molecule)
    img = Draw.MolToImage(mol)

    img.save("molecule.png")
    print("Image saved as molecule.png")

def ADME(molecule):
    mol = Chem.MolFromSmiles(molecule)
    num_rotatable = Lipinski.NumRotatableBonds(mol)

    tpsa_value = Descriptors.TPSA(mol)
    
    print(f"Rotatable bonds: {num_rotatable}")
    print(f"tPSA: {tpsa_value}")

    if tpsa_value <=  140 and num_rotatable <= 10:
        return "Oral Bioavailability status: " + str(True)
    return "Oral Bioavailability status: " + str(False)
