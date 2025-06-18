from rdkit import Chem
from rdkit.Chem import (
    Crippen,
    Descriptors,
    Lipinski,
    PandasTools,
    rdMolDescriptors
)

"""
Script to calculate physicochemical properties of molecules.

The get_physicochemical_properties function calculates 15 molecular descriptors:
- Number of nitrogen atoms
- Number of oxygen atoms  
- Number of chiral centers
- Molecular weight
- Number of heavy atoms
- Number of hydrogen bond acceptors
- Number of hydrogen bond donors
- LogP
- Topological polar surface area
- Number of aromatic atoms
- Number of rings
- Fraction of Csp3 atoms
- Number of sulfur atoms
- Number of halogen atoms
- Molar refractivity

The descriptors are calculated using RDKit's descriptor calculation functions and stored
in separate columns of the input DataFrame.
"""

def get_physicochemical_properties(molDF, smiles_column):
    """
    Applies all property calculations to the ring systems of the dataframe and stores each property in a new column

    Args:
        molDF: dataframe with ring systems as SMILES in the column 'ringSmiles'
        smiles_column: name of column containing SMILES strings
        
    Returns:
        DataFrame with ring system molecules and their properties
    """
    # Convert SMILES to RDKit molecules
    PandasTools.AddMoleculeColumnToFrame(molDF, smiles_column, 'Molecule')
    print('Start calculating properties.')
    
    # Define property calculations grouped by type
    atomic_properties = {
        'N': (get_molecule_composition, 7),
        'O': (get_molecule_composition, 8), 
        'S': (get_molecule_composition, 16)
    }
    
    structural_properties = {
        'chiral': get_nof_chiral_centers,
        'MW': get_MW,
        'heavy_atoms': num_heavy_atoms,
        'numRings': num_rings,
        'frac_csp3': fraction_csp3,
        'numAro': num_aromatic_atoms,
        'nHalogens': num_halogens
    }
    
    physicochemical_properties = {
        'h_acc': (num_of_h_acceptors_and_donors, True),
        'h_don': (num_of_h_acceptors_and_donors, False),
        'logP': get_logp,
        'TPSA': get_TPSA,
        'MR': get_mr
    }
    
    # Calculate atomic properties
    for prop, (func, atomic_num) in atomic_properties.items():
        molDF[prop] = molDF['Molecule'].apply(func, args=(atomic_num,))
        
    # Calculate structural properties  
    for prop, func in structural_properties.items():
        molDF[prop] = molDF['Molecule'].apply(func)
        
    # Calculate physicochemical properties
    for prop, func_args in physicochemical_properties.items():
        if isinstance(func_args, tuple):
            func, arg = func_args
            molDF[prop] = molDF['Molecule'].apply(func, args=(arg,))
        else:
            molDF[prop] = molDF['Molecule'].apply(func_args)


def get_molecule_composition(mol, requestedAtomicNum):
    """
    Counts the number of atoms of a given element in the ring system

    Args:
        mol: the ring system molecule
        requestedAtomicNum: atomic number of the element to count
        
    Returns:
        int: number of atoms of the requested element
    """
    counter = 0
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() == requestedAtomicNum:
            counter += 1
    return counter


def get_nof_chiral_centers(mol):
    """Returns number of chiral centers including unassigned ones"""
    return len(Chem.FindMolChiralCenters(mol, includeUnassigned=True))


def get_MW(mol):
    """Returns molecular weight rounded to 3 decimal places"""
    return round(Descriptors.MolWt(mol), 3)


def num_heavy_atoms(mol):
    """Returns count of non-hydrogen atoms"""
    return Lipinski.HeavyAtomCount(mol)


def num_of_h_acceptors_and_donors(mol, acc=True):
    """Returns count of H-bond acceptors (acc=True) or donors (acc=False)"""
    if acc:
        return Lipinski.NumHAcceptors(mol)
    else:
        return Lipinski.NumHDonors(mol)


def get_logp(mol):
    """Returns calculated LogP value rounded to 3 decimal places"""
    return round(Crippen.MolLogP(mol), 3)


def get_TPSA(mol):
    """Returns topological polar surface area rounded to 3 decimal places"""
    return round(Descriptors.TPSA(mol), 3)


def num_aromatic_atoms(mol):
    """Returns count of aromatic atoms"""
    numAromaticAtoms = 0
    for atom in mol.GetAtoms():
        if atom.GetIsAromatic():
            numAromaticAtoms += 1
    return numAromaticAtoms


def sum_formal_charge(mol):
    """Returns sum of formal charges on all atoms"""
    formalCharge = 0
    for atom in mol.GetAtoms():
        formalCharge += atom.GetFormalCharge()
    return formalCharge


def num_rings(mol):
    """Returns total number of rings"""
    return rdMolDescriptors.CalcNumRings(mol)


def fraction_csp3(mol):
    """Returns fraction of C atoms that are sp3 hybridized"""
    return round(Descriptors.FractionCSP3(mol), 3)


def num_halogens(mol):
    """Returns count of halogen atoms"""
    return Chem.Fragments.fr_halogen(mol)


def get_mr(mol):
    """
    Calculates Wildman-Crippen MR value
    Uses an atom-based scheme based on the values in:
    Wildman and G. M. Crippen JCICS 39 868-873 (1999)
    """
    return round(Crippen.MolMR(mol), 3)
