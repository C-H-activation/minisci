import os

from minisci.helper.features import create_feature_vector
from minisci.helper.protonate import get_dimorphite_protomers
from minisci.helper.rdkit import (
    create_image_from_smiles_barriers,
    extract_aromatic_carbons,
)
from minisci.helper.results import Results
from minisci.helper.scikit import (
    create_boltzmann_weights,
    create_machine_learning_barriers,
)
import typer

app = typer.Typer()


@app.callback()
def main(
    compound: str = typer.Argument(...),
    radical: str = typer.Argument(...),
    number_protons: int = typer.Option(
        0, help="Number of protons to generate protomers."
    ),
    svg_name: str = typer.Option("out.svg", help="Name of the SVG file."),
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
    results = Results()
    results.feature_vectors = dict()
    results.barriers = dict()
    results.weights = dict()
    results.model_path = os.getcwd() + "/minisci/assets/gbr.mod"

    if not number_protons:
        ##################
        # no protomer case
        ##################

        # set aromatic carbon atoms
        results.aromatic_carbons = extract_aromatic_carbons(compound)

        # create feature vectors
        create_feature_vector(results, compound, radical, 0)

        # create machine learning
        results.barriers = create_machine_learning_barriers(
            results.feature_vectors, results.model_path
        )

        # create Boltzmann weights
        results.weights = create_boltzmann_weights(results.barriers)

        # create SVG
        create_image_from_smiles_barriers(
            compound, results.barriers, results.weights, svg_name
        )


if __name__ == "__main__":
    typer.run(main)
