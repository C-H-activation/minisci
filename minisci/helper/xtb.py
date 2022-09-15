"""xTB helper functions."""
import os

from helper.kallisto import kallisto_molecule_from_rdkit_molecule
from helper.rdkit import rdkit_molecule_from_smiles


def smiles_to_xtb_energy(smiles: str) -> float:
    """Obtain GFN2-xtb energy from SMILES string.

    Args:
    smiles: SMILES string for compound (str)

    Returns:
    GFN2-xtb energy in Hartree.
    """
    # construct a molecule
    mol = rdkit_molecule_from_smiles(smiles)
    kmol = kallisto_molecule_from_rdkit_molecule(mol)
    kmol.writeMolecule("kallisto_molecule.xyz")

    # get charge
    charge = int(smiles.count("+"))
    xtb_args = ["--opt", "normal", "--chrg", str(charge), "-P", "1"]
    return perform_xtb_calculation("kallisto_molecule.xyz", xtb_args, os.getcwd())


def perform_xtb_calculation(inp: str, args: list, dirp: str) -> str:
    """Perform an xTB calculation.

    Args:
    inp: String including the input (e.g., structures/input.xyz).
    args: Arguments defining the xTB calculation.
    dirp: Name of temporary directory.
    Returns:
    GFN2-xtb Energy of optimized geoemetry.
    """
    try:
        import subprocess
    except ImportError:
        raise ImportError(
            "workflow requires the module subprocess. Please install the module subprocess."
        )
    try:
        import os
    except ImportError:
        raise ImportError(
            "workflow requires the module os. Please install the module os."
        )

    # save input name
    inputName = extract_input_name(inp)

    # insert elements into arguments
    args.insert(0, "/projects/cp/knkr256/bin/xtb.6.3.2")
    args.insert(1, inputName)

    # get environment settings
    env = os.environ.copy()

    outputName = "xtb.out"
    with open(os.path.join(dirp, outputName), "w", newline=None) as outputFile:
        subprocess.call(
            args,
            shell=False,
            stdin=None,
            stderr=subprocess.STDOUT,
            universal_newlines=False,
            cwd=dirp,
            stdout=outputFile,
            env=env,
        )

    path = os.path.join(dirp, "xtbopt.xyz")
    return extract_xtb_energy(path)


def extract_xtb_energy(inp: str) -> float:
    """Extract energy from optimized xtb structure."""
    if existFile(inp):
        with open(inp, "r") as f:
            lines = f.readlines()
            for i, line in enumerate(lines):
                if i == 1:
                    fields = line.split()
                    energy = fields[1]
                    break
        return float(energy)
    else:
        return float(14)


def extract_input_name(inp: str) -> str:
    """Extract the input name from a path."""

    if os.sep in inp:
        data = inp.split(os.sep)
        inp = data[-1]

    return inp


def existFile(path: str) -> bool:
    """Check if file with 'path' exists."""
    return os.path.isfile(path)
