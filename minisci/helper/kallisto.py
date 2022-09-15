"""kallisto helper functions."""
from kallisto.atom import Atom
from kallisto.molecule import Molecule
from kallisto.units import Bohr
from rdkit import Chem
from rdkit.Chem import GetPeriodicTable


def kallisto_molecule_from_rdkit_molecule(rdkit_molecule: Chem.rdchem.Mol) -> Molecule:
    """Create a kallisto molecule from RDKit molecule.
    Args:
    rdkit_molecule: RDKit molecule
    Returns:
    A kallisto molecule (kallisto.molecule.Molecule)
    Raises:
    KallistoError: An error if the kallisto molecule cannot be created
    """
    # get the name of the molecule if it comes from SDF
    name = ""
    if rdkit_molecule.HasProp("_Name"):
        name = rdkit_molecule.GetProp("_Name")
    # get all xyz coordinates and split into list of lines
    xyz = Chem.rdmolfiles.MolToXYZBlock(rdkit_molecule).split("\n")
    # remove empty lines or molecule name from list
    xyz = [string for string in xyz if string != "" and string != name]
    # remove number of atoms as given in xmol files (first line)
    xyz = xyz[1:]

    # setup periodic table
    pt = GetPeriodicTable()
    # create list of atoms
    atoms = []
    # create kallisto molecule
    for coord in xyz:
        elem, x, y, z = coord.split()[:4]

        # convert atomic coordinates from Angstrom to Bohr
        position = [float(x) / Bohr, float(y) / Bohr, float(z) / Bohr]
        atom = Atom(symbol=pt.GetAtomicNumber(elem), position=position)
        atoms.append(atom)
    kallisto_mol = Molecule(symbols=atoms)
    if "numbers" not in kallisto_mol.arrays.keys():
        raise RuntimeError(
            "The kallisto molecule was not created for the input '{}'".format(
                Chem.MolToSmiles(rdkit_molecule)
            )
        )
    return kallisto_mol
