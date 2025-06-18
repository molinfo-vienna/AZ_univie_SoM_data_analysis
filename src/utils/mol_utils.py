from rdkit import Chem
from rdkit.Chem.Descriptors import descList, MolWt
from rdkit.Chem.Scaffolds import MurckoScaffold
import chembl_structure_pipeline


def get_murcko_scaffold_from_smi(smi):
    return Chem.MolToSmiles(MurckoScaffold.GetScaffoldForMol(Chem.MolFromSmiles(smi)))


def get_murcko_scaffold_from_mol(mol):
    return Chem.MolToSmiles(MurckoScaffold.GetScaffoldForMol(mol))


def get_MW(mol):
    return float(MolWt(mol))


def get_inchi_noStereo_from_smi(smi):
    return Chem.MolToInchi(Chem.MolFromSmiles(Chem.MolToSmiles(Chem.MolFromSmiles(smi),isomericSmiles=False)),options='-SNon')


def get_inchi_noStereo_from_mol(mol):
    return Chem.MolToInchi(Chem.MolFromSmiles(Chem.MolToSmiles(mol,isomericSmiles=False)),options='-SNon')


def get_smi_noStereo_from_mol(mol):
    return Chem.MolToSmiles(mol,isomericSmiles=False)


def get_smi_noStereo_from_smi(smi):
    return Chem.MolToSmiles(Chem.MolFromSmiles(smi),isomericSmiles=False)


def element_label(mol,allowedAtomNrs={1,5,6,7,8,9,14,15,16,17,34,35,53}):
    '''
    allowedAtomNrs={1,5,6,7,8,9,14,15,16,17,34,35,53}
    # H, B, C, N, O, F, Si, S, Cl, Br, Se, I, P
    '''
    
    atomNrs = set()
    for atom in mol.GetAtoms():
        atomNrs.add(atom.GetAtomicNum())
        notAllowedAtomsInMolecule = atomNrs - allowedAtomNrs
    if len(notAllowedAtomsInMolecule) != 0:
        return 0
    else:
        return 1
    
    
def preprocess(mol):
    '''preprocess a molecular structure using chembl_structure_pipline
    '''
    mol.Compute2DCoords()
    mol_new = chembl_structure_pipeline.standardize_mol(mol)
    mol_new,csvErrorCode = chembl_structure_pipeline.get_parent_mol(mol_new)
    
    return mol_new

