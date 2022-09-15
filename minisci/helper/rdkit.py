"""RDKit helper functions."""
import logging
import os

from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import Draw
from rdkit.Chem import Descriptors
from rdkit.Chem.Draw import rdMolDraw2D

from helper.kallisto import kallisto_molecule_from_rdkit_molecule


def rdkit_molecule_from_smiles(smiles: str, minimisation_method=None):
    """Molecule preparation: Parse SMILES, add hydrogens, and does energy minimisation.
    Args:
    smiles: A molecule SMILES string representation (default '')
    Returns:
    An RDKit molecule (rdkit.Chem.rdchem.Mol) or None if the process fails
    """
    # create molecule, add hydrogens, and generate embedding
    m = Chem.MolFromSmiles(smiles)
    if m is None:
        logging.error(
            "The RDKit SMILES parsing has failed for the molecule: %s", smiles
        )
        return None
    mh = Chem.AddHs(m)
    emb_code = AllChem.EmbedMolecule(mh, randomSeed=11)
    if emb_code == -1:
        logging.error("The RDKit embedding has failed for the molecule: %s", smiles)
        return None

    # energy minimisation
    valid_methods = [None, "MMFF94", "MMFF94s", "UFF"]
    if minimisation_method is not None:
        if minimisation_method not in valid_methods:
            raise ValueError(f"Select a valid minimisation method {valid_methods}")
        if minimisation_method == valid_methods[1]:
            AllChem.MMFFOptimizeMolecule(mh)
        if minimisation_method == valid_methods[2]:
            AllChem.MMFFOptimizeMolecule(mh, mmffVariant="MMFF94s")
        if minimisation_method == valid_methods[3]:
            AllChem.UFFOptimizeMolecule(mh)
    return mh


def extract_aromatic_carbons(smiles: str) -> list:
    """Identify all aromatic C-H positions from SMILES.
    Args:
    smiles: SMILES of molecule (str)

    Returns:
    list of aromactic C-H atoms."""

    try:
        from rdkit import Chem
    except ImportError:
        raise ImportError("Module rdkit is required. Please install module rdkit.")
    try:
        import os
    except ImportError:
        raise ImportError("Module os is required. Please install module os.")

    # create substrate from SMILES
    compound = Chem.MolFromSmiles(smiles)

    # identify aromatic Carbon atoms from pattern
    # and unpack tuple to get list of aromatic Carbon atoms
    pattern = Chem.MolFromSmarts("c")
    aromatic_carbons = compound.GetSubstructMatches(pattern)
    aromatic_carbons = [atom for (atom,) in aromatic_carbons]

    # get ring information and extract rings
    ri = compound.GetRingInfo()
    rings = ri.AtomRings()

    # C-H indices
    atoms = []

    for index, atom in enumerate(compound.GetAtoms()):
        for ring in rings:
            #################################################
            # extract Carbon atoms that
            # 1) belong to a ring system
            # 2) are aromatic
            # 3) occur only once (exclude doubles)
            # 4) have ONE C-H bond.
            #################################################
            if (
                (index in ring)
                and (index in aromatic_carbons)
                and (index not in atoms)
                and (atom.GetNumImplicitHs() == 1)
            ):
                atoms.append(index)
    return atoms


def is_smiles_aromatic(smiles: str) -> bool:
    """Check if a SMILES string is aromatic"""
    m = rdkit_molecule_from_smiles(smiles)
    if m is None:
        return False
    ri = m.GetRingInfo()
    for idx, ring in enumerate(ri.BondRings()):
        check = is_ring_aromatic(m, ri.BondRings()[idx])
        # we find an aromatic ring
        if check is True:
            return True
        # there could be many rings -> check next
        continue
    # no aromatic ring found
    return False


def is_ring_aromatic(mol: Chem.rdchem.Mol, bondRing: Chem.rdchem.RingInfo):
    for id in bondRing:
        if not mol.GetBondWithIdx(id).GetIsAromatic():
            return False
    return True


def xmol_file_from_smiles(smiles: str, name: str):
    """Create xmol file with name from SMILES."""
    mol = rdkit_molecule_from_smiles(smiles)
    kmol = kallisto_molecule_from_rdkit_molecule(mol)
    kmol.writeMolecule(f"{name}.xyz")


def extract_radical_idx_from_smiles(radical_smiles: str) -> int:
    """Extract the atomic index of a radical SMILES."""
    rdkit_molecule = rdkit_molecule_from_smiles(radical_smiles)
    # check if not a radical
    if Descriptors.NumRadicalElectrons(rdkit_molecule) == 0:
        return -1
    for idx, atom in enumerate(rdkit_molecule.GetAtoms()):
        number_of_radicals = atom.GetNumRadicalElectrons()
        if number_of_radicals > 0:
            return idx


def create_image_from_smiles_barriers(
    smiles: str,
    barriers: dict,
    weights: dict,
    name: str,
    wdir=os.getcwd(),
):
    """Create a molecular SVG image with RDKit where"""
    # create RDKit molecule
    m = Chem.MolFromSmiles(smiles)
    m2 = Chem.Mol(m)

    # scale barriers: kcal/mol -> kJ/mol
    scaling_factor = 4.184
    barriers_in_kj_mol = {
        key: (val * scaling_factor) for (key, val) in barriers.items()
    }

    # set barriers and weights
    for atom in m2.GetAtoms():
        try:
            idx = atom.GetIdx()
            label_barrier = f"{barriers_in_kj_mol[idx]:0.0f}"
            label_weight = f"{weights[idx]:0.0f}"
            atom.SetProp(
                "atomNote",
                f"{label_barrier} ; {label_weight}",
            )
        except KeyError:
            continue

    # create SVG
    d2d = rdMolDraw2D.MolDraw2DSVG(400, 400)
    d2d.DrawMolecule(m2)
    d2d.FinishDrawing()
    path = os.path.join(wdir, name)
    with open(path, "w") as svg:
        svg.write(d2d.GetDrawingText())
