from typing import Any

from rdkit.Chem import Mol, MolFromMolBlock, MolFromSmiles, MolToMolBlock, MolToSmiles

from rdfreader.exceptions import InvalidMoleculeError
from rdfreader.parse.molblock import get_mol_block_metadata

_UNSET = object()  # sentinel for cached values that have not been computed yet


class Molecule:
    def __init__(
        self,
        mol_block: str = None,
        properties: dict[str, Any] = dict(),
        except_on_invalid_molecule: bool = True,
        component_type: str = None,
    ):
        """Create a molecule object.

        If a mol_block is provided, the molecule will be initialized with that,
        otherwise the molecule will be initialized empty.

        Parameters
        ----------
        mol_block : str
            A mol block string.
        properties : dict[str, Any]
            A dictionary of properties.
        """
        self.except_on_invalid_molecule: bool = except_on_invalid_molecule
        self.properties: dict[str, Any] = properties
        self.mol_block: str = mol_block  # calls mol_block setter
        self.component_type: str = component_type

    @property
    def mol_block(self) -> str:
        """Returns the mol block string of the molecule.

        Returns
        -------
        str
            The mol block string of the molecule.
        """
        return self._mol_block

    @mol_block.setter
    def mol_block(self, mol_block: str):
        """Set the mol block string of the molecule.

        Parameters
        ----------
        mol_block : str
            A mol block string.
        """
        if not mol_block:
            self._mol_block = None
            self._clear_cache()
        else:
            self._from_mol_block(mol_block)

    def _clear_cache(self) -> None:
        """Reset the cached RDKit molecule and SMILES string."""
        self._cached_rd_mol = _UNSET
        self._cached_smiles = _UNSET

    def _get_cached_rd_mol(self) -> Mol:
        """Return the RDKit molecule parsed from the mol block, parsing it only
        once."""
        if self._cached_rd_mol is _UNSET:
            self._cached_rd_mol = MolFromMolBlock(self.mol_block)
        return self._cached_rd_mol

    @property
    def rd_mol(self) -> Mol:
        """Return the RDKit molecule object.

        A new copy is returned on each access, so it can be safely
        modified.
        """
        rd_mol = self._get_cached_rd_mol()
        return Mol(rd_mol) if rd_mol is not None else None

    @property
    def smiles(self) -> str:
        """Returns the SMILES string of the molecule using RDKit.

        Returns:
            str: SMILES string of the molecule.
        """
        if self._cached_smiles is _UNSET:
            try:
                self._cached_smiles = MolToSmiles(self._get_cached_rd_mol())
            except Exception:
                self._cached_smiles = None
        return self._cached_smiles

    @property
    def metadata(self) -> dict[str, Any]:
        """Returns the metadata of the molecule from the mol block string.

        Returns:
            dict[str, Any]: The metadata of the molecule.
        """
        return get_mol_block_metadata(self.mol_block)

    def _from_mol_block(self, mol_block: str, properties: dict[str, Any] = dict()) -> None:
        """Initialize the molecule object with a mol block string.

        Parameters
        ----------
        mol_block : str
            A mol block string.
        """
        self.properties.update(properties)
        self._mol_block = mol_block
        self._clear_cache()

        if self.except_on_invalid_molecule:
            try:
                assert self._get_cached_rd_mol() is not None
            except AssertionError:
                raise InvalidMoleculeError("mol_block is not a valid mol block string.")

    @classmethod
    def from_mol_block(cls, mol_block: str, properties: dict[str, Any] = {}) -> "Molecule":
        """Create a Molecule object from a mol block string.

        Parameters
        ----------
        mol_block : str
            A mol block string.

        Returns
        -------
        Molecule
            A Molecule object.
        """
        mol = cls()
        mol._from_mol_block(mol_block, properties)
        return mol

    def __eq__(self, __o: object) -> bool:
        """Returns True if the molecules are equal.

        Returns
        -------
        bool
            True if the molecules are equal.
        """
        if not isinstance(__o, Molecule):
            return False
        return self.mol_block == __o.mol_block

    def __str__(self) -> str:
        """Returns the smiles string of the molecule.

        Returns
        -------
        str
            The smiles string of the molecule.
        """
        return self.smiles

    def __repr__(self) -> str:
        return f"Molecule(smiles={self.smiles})"

    @classmethod
    def from_smiles(cls, smiles: str) -> "Molecule":
        """Create a Molecule object from a smiles string.

        Parameters
        ----------
        smiles : str
            A smiles string.

        Returns
        -------
        Molecule
            A Molecule object.
        """
        mol_block = MolToMolBlock(MolFromSmiles(smiles))
        mol = cls()
        mol._from_mol_block(mol_block)
        return mol


class Reactant(Molecule):
    def __init__(*args, **kwargs):
        super().__init__(**args, **kwargs, component_type="reactant")


class Product(Molecule):
    def __init__(*args, **kwargs):
        super().__init__(**args, **kwargs, component_type="product")


class Solvent(Molecule):
    def __init__(*args, **kwargs):
        super().__init__(**args, **kwargs, component_type="solvent")


class Catalyst(Molecule):
    def __init__(*args, **kwargs):
        super().__init__(**args, **kwargs, component_type="catalyst")
