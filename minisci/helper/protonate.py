"""Helper function for fast protonation of SMILES strings."""
from typing import List

from dimorphite_dl import DimorphiteDL

from helper.rdkit import (
    is_smiles_aromatic,
)


class Protomer(object):
    """Describe a basic protomer."""

    def __init__(self, smiles: str):
        self.smiles = smiles


def get_dimorphite_protomers(
    input_smiles: str,
    number_protons: int,
    ph: float = 2.0,
) -> List[Protomer]:
    """Get dimorphite protomers at given pH.

    Args:
    smiles: SMILES string representation of molecule (str).
    protons: number of protons to add (int)
    ph: pH value (float).

    Returns:
    List of protomers.
    """
    dimorphite_dl = DimorphiteDL(
        min_ph=ph,
        max_ph=ph + 0.5,
        max_variants=10,
        label_states=False,
        pka_precision=1.0,
    )
    # if we already have "+" in input string
    plus_count = input_smiles.count("+")

    protomers = dimorphite_dl.protonate(input_smiles)

    protonation_sites = [
        protomer
        for protomer in protomers
        if protomer.count("+") == (plus_count + number_protons)
    ]

    protonation_sites = [
        protomer for protomer in protonation_sites if is_smiles_aromatic(protomer)
    ]

    protomers = list()
    for _, smiles in enumerate(protonation_sites):
        if is_smiles_aromatic(smiles):
            protomer = Protomer(smiles=smiles)
            protomers.append(protomer)
    return protomers
