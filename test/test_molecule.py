from unittest.mock import MagicMock

import pytest
from rdkit.Chem import Mol, MolFromMolBlock, MolFromSmiles, MolToMolBlock, MolToSmiles

from rdfreader.chem.mol import Molecule


def assert_molecule_from_mol_block(mol: Molecule, mol_block: str):
    """Assert that a molecule object created from a mol block string has the
    correct properties."""
    assert mol.rd_mol is not None
    assert mol.mol_block == mol_block


def test_molecule_from_mol_block(sample_mol_block):
    """Test that the Molecule.from_mol_block method works, creating a Molecule
    object from a mol block string and correctly creates class attributes."""
    mol = Molecule.from_mol_block(sample_mol_block)
    assert_molecule_from_mol_block(mol, sample_mol_block)


def test_molecule_init_from_mol_block(sample_mol_block):
    """Test we can init a molecule object directly by passing a mol block
    string."""
    mol = Molecule(sample_mol_block)
    assert_molecule_from_mol_block(mol, sample_mol_block)


def test_create_empty_molecule():
    """Test that a Molecule object can be created without a mol block
    string."""
    mol = Molecule()
    assert mol.mol_block is None


def test_molecule_to_rdkit_mol(sample_molecule, sample_mol_block):
    """Test the Molecule.rd_mol property."""
    rd_mol: Mol = sample_molecule.rd_mol
    # can't directly compare Mol objects, so we'll just check that it is of
    # the right type and is not None
    assert rd_mol is not None
    assert isinstance(rd_mol, Mol)


def test_molecule_to_smiles(sample_molecule, sample_mol_block):
    """Test the output of the Molecule.smiles property.

    Verifies the smiles matches that output by RDKit.
    """
    rd_mol: Mol = MolFromMolBlock(sample_mol_block)
    rd_smiles: str = MolToSmiles(rd_mol)
    assert sample_molecule.smiles == rd_smiles


def test_molecule_metadata(sample_molecule, sample_molecule_metadata):
    """Test the Molecule.metadata property."""
    assert sample_molecule.metadata == sample_molecule_metadata


def test_molecule_from_smiles():
    """Test the Molecule.from_smiles method."""
    smiles = "OCO"
    mol = Molecule.from_smiles(smiles)
    assert mol.rd_mol is not None
    assert mol.smiles == smiles


def test_molecule_mol_block_parsed_once(mocker, sample_mol_block):
    """Test that validation and smiles share a single RDKit parse of the mol
    block."""
    parse_spy: MagicMock = mocker.patch("rdfreader.chem.mol.MolFromMolBlock", wraps=MolFromMolBlock)
    mol = Molecule(sample_mol_block)  # validates the mol block
    mol.smiles
    mol.smiles
    parse_spy.assert_called_once_with(sample_mol_block)


def test_molecule_rd_mol_not_cached(mocker, sample_mol_block):
    """Test that the RDKit molecule is not kept in memory: each rd_mol access
    re-parses the mol block."""
    parse_spy: MagicMock = mocker.patch("rdfreader.chem.mol.MolFromMolBlock", wraps=MolFromMolBlock)
    mol = Molecule(sample_mol_block)
    mol.rd_mol
    mol.rd_mol
    assert parse_spy.call_count == 3  # validation + 2 accesses
    assert "_smiles" in mol.__dict__
    assert not any(isinstance(value, Mol) for value in mol.__dict__.values())


def test_molecule_smiles_computed_once(mocker, sample_mol_block):
    """Test that the SMILES string is only generated once by RDKit."""
    smiles_spy: MagicMock = mocker.patch("rdfreader.chem.mol.MolToSmiles", wraps=MolToSmiles)
    mol = Molecule(sample_mol_block)
    first_smiles = mol.smiles
    assert mol.smiles == first_smiles
    smiles_spy.assert_called_once()


def test_molecule_invalid_mol_block_parsed_once(mocker):
    """Test that a failed SMILES computation is cached too, rather than the
    parse being retried on each smiles access (rd_mol always re-parses)."""
    parse_spy: MagicMock = mocker.patch("rdfreader.chem.mol.MolFromMolBlock", wraps=MolFromMolBlock)
    mol = Molecule("not a mol block", except_on_invalid_molecule=False)
    assert mol.smiles is None
    assert mol.smiles is None
    parse_spy.assert_called_once()
    assert mol.rd_mol is None


def test_molecule_rd_mol_is_independent(sample_molecule):
    """Test that modifying the returned RDKit molecule does not affect the
    Molecule, its later rd_mol values or its cached SMILES string."""
    original_smiles = sample_molecule.smiles
    rd_mol = sample_molecule.rd_mol
    assert rd_mol is not sample_molecule.rd_mol

    rd_mol.GetAtomWithIdx(0).SetAtomicNum(54)
    assert MolToSmiles(rd_mol) != original_smiles
    assert MolToSmiles(sample_molecule.rd_mol) == original_smiles
    assert sample_molecule.smiles == original_smiles


@pytest.mark.parametrize(
    "new_mol_block, expected_smiles",
    [
        (MolToMolBlock(MolFromSmiles("OCO")), "OCO"),
        (None, None),
        ("", None),
    ],
)
def test_molecule_cache_cleared_on_mol_block_change(sample_molecule, new_mol_block, expected_smiles):
    """Test that setting a new mol block, or clearing it, discards the cached
    SMILES string."""
    assert sample_molecule.smiles is not None  # populate the cache
    sample_molecule.mol_block = new_mol_block
    assert sample_molecule.smiles == expected_smiles
    if expected_smiles is not None:
        assert MolToSmiles(sample_molecule.rd_mol) == expected_smiles
