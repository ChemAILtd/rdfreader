import logging
from typing import Any, Optional

from rdkit.Chem.rdChemReactions import ChemicalReaction, ReactionFromSmarts

from rdfreader.chem.mol import Molecule
from rdfreader.chem.utils import reaction_smiles
from rdfreader.exceptions import InvalidReactionError
from rdfreader.parse.rxnblock import DatumParser, get_rxn_block_metadata, mol_blocks_from_rxn_block, validate_rxn_block
from rdfreader.parse.utils import CTF_RXNBLOCK_HEADER_FORMAT_STRING

logger = logging.getLogger(__name__)


class Reaction:
    def __init__(
        self,
        rxn_block: str = None,
        id: str = None,
        rdf_metadata: Optional[dict[str, Any]] = None,  # metadata from the rdf file the block was contained in,
        except_on_invalid_molecule: bool = True,
        header_format_string: str = CTF_RXNBLOCK_HEADER_FORMAT_STRING,
        lineno: Optional[int] = None,
        rdf_file: Optional[str] = None,
    ):
        """Create a reaction object.

        If a reaction block is provided, the reaction will be
        initialized with that, otherwise the reaction will be
        initialized empty.
        """
        self.rxn_block = rxn_block
        self.rdf_metadata = rdf_metadata
        self.id = id
        self.lineno = lineno
        self.rdf_file = rdf_file

        self.products: list[Molecule] = list()
        self.reactants: list[Molecule] = list()
        self.catalysts: list[Molecule] = list()
        self.solvents: list[Molecule] = list()
        self.other_reagents: list[Molecule] = list()

        self.properties: dict[str, str] = dict()
        self.metadata: dict[str, Any] = dict()

        if self.rxn_block is not None:
            self._from_rxn_block(
                header_format_string=header_format_string,
                except_on_invalid_molecule=except_on_invalid_molecule,
            )

    def __str__(self) -> str:
        return f"Reaction({self.id}, {self.smiles})"

    def _from_rxn_block(
        self,
        header_format_string: str = CTF_RXNBLOCK_HEADER_FORMAT_STRING,
        except_on_invalid_molecule: bool = True,
    ) -> None:
        """Initialize the reaction object from a reaction block.

        Parameters
        ----------
        rdf_path : Path
            A path to a rdf file.
        """
        if not validate_rxn_block(self.rxn_block):
            raise ValueError("Reaction block is invalid.")

        # get the headers from the rxn block.
        self.metadata.update(get_rxn_block_metadata(self.rxn_block, header_format_string=header_format_string))

        # get the reactants and product mol_blocks, create a molecule object
        # for each, and add it to the reaction.
        reactant_mol_blocks, product_mol_blocks = mol_blocks_from_rxn_block(
            self.rxn_block,
            self.metadata["reactant_count"],
            self.metadata["product_count"],
        )
        for mol_block in reactant_mol_blocks:
            self.reactants.append(
                Molecule(
                    mol_block,
                    except_on_invalid_molecule=except_on_invalid_molecule,
                )
            )

        for mol_block in product_mol_blocks:
            self.products.append(
                Molecule(
                    mol_block,
                    except_on_invalid_molecule=except_on_invalid_molecule,
                )
            )

        # use functions in rdfreader/parse/rdf.py to pull out dtype/datum
        # pairs from the rxn block
        datum_parser = DatumParser(
            self.rxn_block,
        )

        for dtype, datum in datum_parser:
            if not isinstance(datum, Molecule):
                self.properties[dtype] = datum
            else:
                # add to the appropriate molecule list
                getattr(self, f"{datum.component_type}s").append(datum)

        if self.rd_rxn is None:
            raise ValueError("Invalid reaction: couldn't parse in rdkit.")

    @property
    def smiles(self) -> str:
        return reaction_smiles(
            self.reactants,
            self.products,
            self.reagents,
        )

    @property
    def smiles_no_reagents(self) -> str:
        return reaction_smiles(
            self.reactants,
            self.products,
        )

    @property
    def rd_rxn(self) -> ChemicalReaction:
        try:
            return ReactionFromSmarts(self.smiles, useSmiles=True)
        except ValueError as e:
            raise InvalidReactionError(f"Invalid reaction: {e}") from e

    @property
    def reagents(self) -> list[Molecule]:
        """Return a single list of all reagents.

        Returns
        -------
        list[Molecule]
            A list of all reagents.
        """
        return self.catalysts + self.solvents + self.other_reagents
