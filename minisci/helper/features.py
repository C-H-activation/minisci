"""Featurization helper functions."""

from morfeus import BuriedVolume
from morfeus import read_xyz
from morfeus import XTB

try:
    from openbabel import pybel
except ImportError:
    raise ImportError(
        "Please install openbabel using conda (run `conda install -c conda-forge openbabel`)."
    )

from helper.bash import remove_file
from helper.graph import Graph
from helper.kallisto import kallisto_molecule_from_rdkit_molecule
from helper.rdkit import (
    extract_radical_idx_from_smiles,
    rdkit_molecule_from_smiles,
    xmol_file_from_smiles,
)
from helper.results import Results


def create_feature_vector(results: Results, compound: str, radical: str):
    """Create physicochemical feature vector for compound, radical, and barrier."""
    # print information and dump compound
    print("Featurization applied for")
    print(f"Compound: {compound}")
    print(f"Radical: {radical}")
    m = rdkit_molecule_from_smiles(compound)
    km = kallisto_molecule_from_rdkit_molecule(m)
    km.writeMolecule("compound_features.xyz")

    # get aromatic carbon sites
    sites = results.aromatic_carbons

    # initialize features
    features = {}

    # set SMILES
    compound_smiles = compound
    radical_smiles = radical

    # molecular charge
    molecular_charge = 0

    # extract charge of input structure
    plus_count = compound_smiles.count("+")
    minus_count = compound_smiles.count("-")
    if plus_count > 0:
        molecular_charge = plus_count - minus_count

    # compound
    xmol_file_from_smiles(compound_smiles, "compound")
    compound_elements, compound_coordinates = read_xyz("compound.xyz")
    compound_xtb = XTB(
        elements=compound_elements,
        coordinates=compound_coordinates,
        charge=molecular_charge,
    )
    compound_xtb_cation = XTB(
        elements=compound_elements,
        coordinates=compound_coordinates,
        charge=molecular_charge + 1,
    )
    compound_xtb_anion = XTB(
        elements=compound_elements,
        coordinates=compound_coordinates,
        charge=molecular_charge - 1,
    )

    # compound partial charges
    compound_charges = compound_xtb.get_charges()
    compound_charges_cation = compound_xtb_cation.get_charges()
    compound_charges_anion = compound_xtb_anion.get_charges()
    results.charges = compound_charges

    # compound Fukui coefficients
    compound_electrophilicity = compound_xtb.get_fukui("electrophilicity")
    compound_local_electrophilicity = compound_xtb.get_fukui("local_electrophilicity")
    compound_nucleophilicity = compound_xtb.get_fukui("nucleophilicity")
    compound_local_nucleophilicity = compound_xtb.get_fukui("local_nucleophilicity")
    results.electrophilicity = compound_electrophilicity
    results.local_electrophilicity = compound_local_electrophilicity
    results.nucleophilicity = compound_nucleophilicity
    results.local_nucleophilicity = compound_local_nucleophilicity

    # other electronical features of compound
    compound_ip = compound_xtb.get_ip(corrected=True)
    compound_ea = compound_xtb.get_ea()
    compound_mu = (compound_ea - compound_ip) / 2.0
    compound_eta = compound_ip - compound_ea
    compound_homo = compound_xtb.get_homo()
    compound_lumo = compound_xtb.get_lumo()

    # radical
    radical_idx = extract_radical_idx_from_smiles(radical_smiles)

    # sort radical if necessary
    if radical_idx != 0:
        radical_mol = rdkit_molecule_from_smiles(radical_smiles)
        radical_kmol = kallisto_molecule_from_rdkit_molecule(radical_mol)
        radical_bonds = radical_kmol.get_bonds()
        radical_nat = radical_kmol.get_number_of_atoms()

        # create molecular graph from kallisto molecule
        g = Graph(radical_kmol)
        for i in range(radical_nat):
            partners = radical_bonds[i]
            for j in partners:
                # add edges
                g.addEdge(j, i)
        g.BFS(int(radical_idx), f"radical.xyz")
        # get new radical mol from sorted xyz
        mol = next(pybel.readfile("xyz", f"radical.xyz"))
        radical_smiles = mol.write(format="smi")
        radical_mol = rdkit_molecule_from_smiles(radical_smiles)
        radical_kmol = kallisto_molecule_from_rdkit_molecule(radical_mol)
    else:
        radical_mol = rdkit_molecule_from_smiles(radical_smiles)
        radical_kmol = kallisto_molecule_from_rdkit_molecule(radical_mol)
        radical_kmol.writeMolecule("radical.xyz")

    radical_elements, radical_coordinates = read_xyz("radical.xyz")
    radical_xtb = XTB(
        elements=radical_elements, coordinates=radical_coordinates, n_unpaired=1
    )
    radical_xtb_cation = XTB(
        elements=radical_elements,
        coordinates=radical_coordinates,
        n_unpaired=1,
        charge=1,
    )
    radical_xtb_anion = XTB(
        elements=radical_elements,
        coordinates=radical_coordinates,
        n_unpaired=1,
        charge=-1,
    )

    # radical partial charges
    radical_charges = radical_xtb.get_charges()
    radical_charges_cation = radical_xtb_cation.get_charges()
    radical_charges_anion = radical_xtb_anion.get_charges()

    # radical Fukui coefficients
    radical_electrophilicity = radical_xtb.get_fukui("electrophilicity")
    radical_local_electrophilicity = radical_xtb.get_fukui("local_electrophilicity")
    radical_nucleophilicity = radical_xtb.get_fukui("nucleophilicity")
    radical_local_nucleophilicity = radical_xtb.get_fukui("local_nucleophilicity")

    # other electronical features of radical
    radical_ip = radical_xtb.get_ip(corrected=True)
    radical_ea = radical_xtb.get_ea()
    radical_mu = (radical_ea - radical_ip) / 2.0
    radical_eta = radical_ip - radical_ea
    radical_homo = radical_xtb.get_homo()
    radical_lumo = radical_xtb.get_lumo()

    # loop over all aromatic carbon atoms
    for _, carbon_idx in enumerate(sites):
        # we set the reactive site instead of obtaining it from the product
        reactive_site = carbon_idx

        # Morfeus returns indices starting by 1 instead of 0
        reactive_site += 1
        # radical geometry is sorted and radical atom comes first
        radical_site = 1

        # compound buried volume
        compound_bv3 = BuriedVolume(
            compound_elements, compound_coordinates, reactive_site, radius=3
        )
        compound_bv4 = BuriedVolume(
            compound_elements, compound_coordinates, reactive_site, radius=4
        )
        compound_bv5 = BuriedVolume(
            compound_elements, compound_coordinates, reactive_site, radius=5
        )
        compound_bv3_fraction = compound_bv3.fraction_buried_volume
        compound_bv4_fraction = compound_bv4.fraction_buried_volume
        compound_bv5_fraction = compound_bv5.fraction_buried_volume

        # radical buried volume
        radical_bv3 = BuriedVolume(
            radical_elements, radical_coordinates, radical_site, radius=3
        )
        radical_bv4 = BuriedVolume(
            radical_elements, radical_coordinates, radical_site, radius=4
        )
        radical_bv5 = BuriedVolume(
            radical_elements, radical_coordinates, radical_site, radius=5
        )
        radical_bv3_fraction = radical_bv3.fraction_buried_volume
        radical_bv4_fraction = radical_bv4.fraction_buried_volume
        radical_bv5_fraction = radical_bv5.fraction_buried_volume

        feature_vector = []
        # start with compound
        # 1. charge carbon reactive site neutral
        feature_vector.append(compound_charges[reactive_site])
        # 2. charge carbon reactive site cation
        feature_vector.append(compound_charges_cation[reactive_site])
        # 3. charge carbon reactive site anion
        feature_vector.append(compound_charges_anion[reactive_site])
        # 4. ionization potential
        feature_vector.append(compound_ip)
        # 5. electron affinity
        feature_vector.append(compound_ea)
        # 6. HOMO energy
        feature_vector.append(compound_homo)
        # 7. LUMO energy
        feature_vector.append(compound_lumo)
        # 8. carbon reactive site electrophilicity (Fukui)
        feature_vector.append(compound_electrophilicity[reactive_site])
        # 9. carbon reactive site nucleophilicity (Fukui)
        feature_vector.append(compound_nucleophilicity[reactive_site])
        # 10. carbon reactive site local electrophilicity (Fukui)
        feature_vector.append(compound_local_electrophilicity[reactive_site])
        # 11. carbon reactive site local nucleophilicity (Fukui)
        feature_vector.append(compound_local_nucleophilicity[reactive_site])
        # 12. chemical potential
        feature_vector.append(compound_mu)
        # 13. chemical hardness
        feature_vector.append(compound_eta)
        # 14. fraction buried volume radius = 3 Angstrom
        feature_vector.append(compound_bv3_fraction)
        # 15. fraction buried volume radius = 4 Angstrom
        feature_vector.append(compound_bv4_fraction)
        # 16. fraction buried volume radius = 5 Angstrom
        feature_vector.append(compound_bv5_fraction)

        # continue with radical
        # 17. charge carbon radical
        feature_vector.append(radical_charges[radical_site])
        # 18. charge carbon radical cation
        feature_vector.append(radical_charges_cation[radical_site])
        # 19. charge carbon radical anion
        feature_vector.append(radical_charges_anion[radical_site])
        # 20. ionization potential
        feature_vector.append(radical_ip)
        # 21. electron affinity
        feature_vector.append(radical_ea)
        # 22. HOMO energy
        feature_vector.append(radical_homo)
        # 23. LUMO energy
        feature_vector.append(radical_lumo)
        # 24. carbon radical electrophilicity (Fukui)
        feature_vector.append(radical_electrophilicity[radical_site])
        # 25. carbon radical nucleophilicity (Fukui)
        feature_vector.append(radical_nucleophilicity[radical_site])
        # 26. carbon radical local electrophilicity (Fukui)
        feature_vector.append(radical_local_electrophilicity[radical_site])
        # 27. carbon radical local nucleophilicity (Fukui)
        feature_vector.append(radical_local_nucleophilicity[radical_site])
        # 28. chemical potential
        feature_vector.append(radical_mu)
        # 29. chemical hardness
        feature_vector.append(radical_eta)
        # 30. fraction buried volume radius = 3 Angstrom
        feature_vector.append(radical_bv3_fraction)
        # 31. fraction buried volume radius = 4 Angstrom
        feature_vector.append(radical_bv4_fraction)
        # 32. fraction buried volume radius = 5 Angstrom
        feature_vector.append(radical_bv5_fraction)
        features[carbon_idx] = feature_vector

    # clean up directory
    remove_file(f"compound.xyz")
    remove_file(f"radical.xyz")
    remove_file("charges")
    remove_file("xtbtopo.mol")
    remove_file("xtbrestart")
    remove_file("wbo")
    remove_file("charges")
    results.feature_vectors = features
