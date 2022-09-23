import os

import typer

from minisci.helper.features import create_feature_vector
from minisci.helper.protonate import get_dimorphite_protomers
from minisci.helper.rdkit import (
    create_image_from_smiles_barriers,
    extract_aromatic_carbons,
    standardize,
)
from minisci.helper.results import Results
from minisci.helper.scikit import (
    create_boltzmann_weights,
    create_machine_learning_barriers,
)

app = typer.Typer()


@app.callback()
def main(
    compound: str = typer.Argument(...),
    radical: str = typer.Argument(...),
    number_protons: int = typer.Option(
        0, help="Number of protons to generate protomers."
    ),
    svg_name: str = typer.Option("compound_image", help="Name of the SVG file."),
):
    """
    Regioselectivity determination for Minisci reactions.

    Args:
    compound: compound SMILES string (str)
    radical: radical SMILES string (str)
    protons: number of protons (int, default=0)

    """

    # typer.echo(f"Run regioselectivity determination for {compound}")
    # typer.echo(f"Run regioselectivity determination with radical {radical}")
    # typer.echo(f"Run regioselectivity determination with protons {number_protons}")

    # generate protomers if necessary
    protomers = []
    if number_protons:
        protomers = get_dimorphite_protomers(compound, number_protons)

    # initialize results object
    results = initialise_results_object()

    if not number_protons:
        ##################
        # no protomer case
        ##################

        # standardize SMILES
        compound = standardize(compound)

        # prepare the results object
        results = prepare_results_object_for_compound(results, compound, radical)
        print("Compound, radical:", compound, radical)

        # create SVG for compound
        create_image_from_smiles_barriers(
            compound, results.barriers, results.weights, f"{svg_name}.svg"
        )
    else:
        for idx, protomer in enumerate(protomers):
            # set SVG name
            svg_name = "protomer_image"

            # extract protomer compound
            protomer_compound = protomer.smiles

            # standardize protomer compound
            protomer_compound = standardize(protomer_compound)

            # initialize results object
            results = initialise_results_object()

            print("Compound, radical:", protomer_compound, radical)
            # prepare results object
            results = prepare_results_object_for_compound(
                results, protomer_compound, radical
            )

            # create SVG for protomer
            create_image_from_smiles_barriers(
                protomer_compound,
                results.barriers,
                results.weights,
                f"{svg_name}_{idx}.svg",
            )


def initialise_results_object() -> Results:
    """Initialises an empty results object."""
    results = Results()
    results.feature_vectors = dict()
    results.barriers = dict()
    results.weights = dict()
    results.model_path = os.getcwd() + "/minisci/assets/gbr.mod"
    return results


def prepare_results_object_for_compound(
    results: Results, compound: str, radical: str
) -> Results:
    """Prepare results object for a given compound and radical SMILES."""
    # set aromatic carbon atoms in results object
    results.aromatic_carbons = extract_aromatic_carbons(compound)

    # create feature vectors
    create_feature_vector(results, compound, radical)

    # create machine learning
    results.barriers = create_machine_learning_barriers(
        results.feature_vectors, results.model_path
    )

    # create Boltzmann weights
    results.weights = create_boltzmann_weights(results.barriers)

    return results


if __name__ == "__main__":
    typer.run(main)
