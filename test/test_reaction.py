# Typing
from unittest.mock import MagicMock

import pytest
from rdkit.Chem import MolFromMolBlock

from rdfreader.chem.reaction import Reaction


def test_reaction_from_rxn_block(sample_rxn_block, sample_rxn_block_metadata):
    reaction = Reaction(sample_rxn_block)

    assert isinstance(reaction, Reaction)
    assert reaction.metadata == sample_rxn_block_metadata
    assert reaction.rxn_block == sample_rxn_block
    assert reaction.rd_rxn is not None
    # TODO: check reaction properties match


def test_reaction_from_rxn_block_invalid_raises():
    with pytest.raises(ValueError):
        Reaction("invalid rxn block")


def test_reaction_from_rxn_block_empty_raises():
    with pytest.raises(ValueError):
        Reaction("")


def test_reaction_to_smiles(mocker, sample_rxn_block):
    """Test that the reaction_to_smiles function gets called."""
    reaction_smiles_patch: MagicMock = mocker.patch("rdfreader.chem.reaction.reaction_smiles", return_value="CC>>CC")
    # give a dummy smiles just so the validation check within
    # Reaction.from_rxn_block() passes
    reaction: Reaction = Reaction(sample_rxn_block)
    reaction.smiles
    reaction_smiles_patch.assert_called()
    # not _once because the reaction_to_smiles function is also called when
    # the Reaction object is instantiated


def test_reaction_to_smiles_no_reagents(mocker, sample_rxn_block):
    """Test that the reaction_to_smiles function gets called."""
    reaction_smiles_patch: MagicMock = mocker.patch("rdfreader.chem.reaction.reaction_smiles", return_value="CC>>CC")
    # give a dummy smiles just so the validation check within
    # Reaction.from_rxn_block() passes
    reaction: Reaction = Reaction(sample_rxn_block)
    reaction.smiles_no_reagents
    reaction_smiles_patch.assert_called()


def test_reaction_parses_each_mol_block_once(mocker, sample_rxn_block):
    """Test that parsing a reaction (including validation and building the
    reaction SMILES) only parses each molecule's mol block once with RDKit."""
    parse_spy: MagicMock = mocker.patch("rdfreader.chem.mol.MolFromMolBlock", wraps=MolFromMolBlock)
    reaction: Reaction = Reaction(sample_rxn_block)
    reaction.smiles
    reaction.smiles_no_reagents
    molecule_count = len(reaction.reactants + reaction.products + reaction.reagents)
    assert molecule_count > 0
    assert parse_spy.call_count == molecule_count
