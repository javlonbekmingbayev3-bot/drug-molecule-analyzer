from rdkit.Chem import Descriptors, Lipinski, Draw
from rdkit import Chem

def molecular_weight(molecule):
    """Calculates the molecular weight of the given molecule in SMILES format."""
    m_w = Descriptors.MolWt(molecule)
    return m_w

def calculate_solubility(molecule):
    """Calculates the solubility of the molecule."""
    logP = Descriptors.MolLogP(molecule)
    if logP<1:
        return 'very soluble'
    elif logP> 1 and logP<3:
        return 'moderate soluble'
    else:
        return 'poorly soluble'
    
def ro_5(molecule):
    """Applies Lipinski's rule of five to test whether the molecule is orally bioavailable to the human body."""
    hbd_count = Lipinski.NumHDonors(molecule)
    hba_count = Lipinski.NumHAcceptors(molecule)
    m_w = molecular_weight(molecule)
    LogP = Descriptors.MolLogP(molecule)

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
    img = Draw.MolToImage(molecule)

    img.save("molecule.png")
    print("Image saved as molecule.png\n")

def ADME(molecule):
    """Calculates number of rotatable bonds and Total Polar Surface Area of the molecule."""
    num_rotatable = Lipinski.NumRotatableBonds(molecule)

    tpsa_value = Descriptors.TPSA(molecule)
    
    print(f"Rotatable bonds: {num_rotatable}")
    print(f"tPSA: {tpsa_value}")

    if tpsa_value <=  140 and num_rotatable <= 10:
        return "Oral Bioavailability status: " + str(True)
    return "Oral Bioavailability status: " + str(False)
