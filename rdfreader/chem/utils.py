from rdfreader.chem.mol import Molecule


def mol_list_to_smiles(mol_list: list[Molecule]) -> str:
    """Convert a list of molecules to a SMILES string.

    Parameters
    ----------
    mol_list : list[Molecule]
        The list of molecules to convert to a SMILES string.

    Returns
    -------
    str
        The SMILES string.
    """
    return ".".join([mol.smiles for mol in mol_list if mol.smiles is not None])


def reaction_smiles(
    reactants: list[Molecule],
    products: list[Molecule],
    reagents: list[Molecule] = [],
) -> str:
    """Create a reaction smiles string from lists of product, reactant, and
    reagent molecules."""
    product_smiles = mol_list_to_smiles(products)
    reactant_smiles = mol_list_to_smiles(reactants)
    reagent_smiles = mol_list_to_smiles(reagents)

    return ">".join([reactant_smiles, reagent_smiles, product_smiles])
