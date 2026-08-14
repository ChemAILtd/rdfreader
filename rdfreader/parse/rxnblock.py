from __future__ import annotations

import logging
import re
from typing import Any, Callable, Iterator

from rdfreader.chem.mol import Molecule
from rdfreader.parse.utils import (
    CTF_COMPONENT_COUNT_FORMAT_STRING,
    CTF_DEFAULT_LETTER_TO_FIELD_MAPPING,
    CTF_RXNBLOCK_HEADER_FORMAT_STRING,
    _parse_block_header_line,
    get_whole_line_item,
    make_string_python_safe,
)

logger = logging.getLogger(__name__)


def get_rxn_block_metadata(
    rxn_block: str,
    header_format_string: str = CTF_RXNBLOCK_HEADER_FORMAT_STRING,
    header_field_mapping: dict[str, tuple[str, Callable, Any]] = CTF_DEFAULT_LETTER_TO_FIELD_MAPPING,
    reactant_product_count_format_string: str = CTF_COMPONENT_COUNT_FORMAT_STRING,  # noqa F401
) -> dict[str, Any]:
    """Get the metadata from a reaction block.

    Parameters
    ----------
    rxn_block : str
        A reaction block string.

    Returns
    -------
    dict
        A dictionary of metadata.
    """
    metadata = {}
    rxn_block_lines: list[str] = rxn_block.split("\n")
    metadata["reaction_name"] = get_whole_line_item(rxn_block_lines[1])
    metadata.update(_parse_block_header_line(rxn_block_lines[2], header_format_string, header_field_mapping))
    metadata["comment"] = get_whole_line_item(rxn_block_lines[3])
    metadata.update(
        _parse_block_header_line(
            rxn_block_lines[4],
            reactant_product_count_format_string,
            header_field_mapping,
        )
    )

    return metadata


class DatumParser:
    """Process dtype/datum pairs.

    The class processes the dtype string into a callable method name and then
    calls that method with the datum string to process the data.

    To use: instantiate the class and then call the class as a function with
    the dtype and datum strings.

    Example:
    >>> dp = DatumParser()
    >>> dp("$DATUM", "RXNVARIATION:STEPNO:SOLVENT:MOLSYMBOL")
    """

    def __init__(self, rxn_block: str, except_on_invalid_molecule: bool = True):
        self.rxn_block = rxn_block
        self.except_on_invalid_molecule = except_on_invalid_molecule

    def __iter__(self) -> Iterator[tuple[Any, int, str]]:
        """Iterate over the lines in the reaction block, identifying
        dtype/datum pairs and returning them as a tuple."""
        dtype_string_identifier = "$DTYPE "

        lines = self.rxn_block.split("\n")

        for line_idx, line in enumerate(lines):
            if line.startswith(dtype_string_identifier):
                dtype_string = line
                datum_string = ""
                # capture all the following lines until we find a line that
                # starts with $DTYPE again
                for datum_line in lines[line_idx + 1 :]:
                    if datum_line.startswith(dtype_string_identifier):
                        break

                    if datum_line.endswith("+"):
                        # if the line ends with a plus sign, it is a
                        # continuation of the previous line
                        datum_string += datum_line[:-1]
                    else:
                        datum_string += datum_line + "\n"

                yield self.parse_datum(dtype_string, datum_string)

    def __call__(self, *args, **kwargs):
        """Wraps parse_datum method."""
        return self.parse_datum(*args, **kwargs)

    def parse_datum(self, dtype: str, datum: str) -> tuple[str, str | Molecule]:
        """Process the $DATUM datum.

        Parameters
        ----------
        dtype : str
            The $DTYPE string from the reaction block.
        datum : str
            The corresponding datum string from the reaction block.

        Returns
        -------
        tuple[Any, str]
            A tuple of the processed datum and the parsed
            dtype string.
        """
        parsed_dtype = parse_dtype_string(dtype)
        datum = preprocess_datum_string(datum)

        # if the datum_string is a molblock:
        if detect_molblock_from_datum(datum):
            datum = preprocess_datum_molblock(datum)
            # try and infer the type from the dtype string, etc.)
            reagent_type = "other_reagent"
            for _reagent_type in ["catalyst", "solvent"]:
                if _reagent_type in parsed_dtype or _reagent_type.upper() in parsed_dtype:
                    reagent_type = _reagent_type

            datum = Molecule(
                mol_block=datum, component_type=reagent_type, except_on_invalid_molecule=self.except_on_invalid_molecule
            )
        else:
            datum = datum.strip()

        return parsed_dtype, datum


def detect_molblock_from_datum(datum: str) -> bool:
    """Tries to detect whether a datum string is a molblock."""
    if datum.startswith("$MFMT"):
        return True
    return False


def preprocess_datum_molblock(datum: str) -> str:
    """Strips the datum down to just the molblock."""
    return "\n".join(datum.splitlines()[1:])  # remove the $MFMT line


def preprocess_datum_string(datum: str) -> str:
    """Preprocess the datum string.

    Parameters
    ----------
    datum_str : str
        The datum string.

    Returns
    -------
    str
        The preprocessed datum string.
    """
    re_pattern = r"^\$DATUM\s?"
    datum = re.sub(re_pattern, "", datum)
    return datum


def parse_dtype_string(dtype_string: str) -> str:
    """Parse a $dtype line from a reaction block.

    Returns the contents of the dtype line, with the prefix "$DTYPE " removed.

    Parameters
    ----------
    dtype_string : str
        A string containing a $dtype line.

    Returns
    -------
    str
        The contents of the $dtype line.
    """
    dtype_string = dtype_string.strip()  # Remove leading and trailing whitespace.
    dtype_string = dtype_string.replace("$DTYPE ", "")  # Remove the $dtype tag.
    dtype_string = make_string_python_safe(dtype_string)
    return dtype_string


def validate_rxn_block(rxn_block: str) -> bool:
    """Validates a rxn block.

    Parameters
    ----------
    rxn_block : str
        The rxn block.

    Returns
    -------
    bool
        True if the rxn block is valid.
    """
    if rxn_block.startswith("$RXN"):
        return True
    else:
        return False


def mol_blocks_from_rxn_block(rxn_block: str, reactant_count: int, product_count: int) -> tuple[list[str], list[str]]:
    """Get the mol blocks corresponding to reactants and products from the rxn
    block.

    Params
    ------
    rxn_block : str
        A reaction block string.
    reactant_count : int
        The number of reactants in the reaction.
    product_count : int
        The number of products in the reaction.

    Returns
    -------
    reactants : list[str]
        The mol blocks corresponding to reactants.
    products : list[str]
        The mol blocks corresponding to products.
    """
    reactants: list[str] = []
    products: list[str] = []

    start_string = "$MOL\n"
    end_string = "M  END\n"

    # mol_blocks start with "$MOL" and end with "M  END"
    # so we can just split on these and get the reactants and products.
    # reactants are first, then products.
    mol_blocks: list[str] = rxn_block.split(start_string)

    # remove the first element as it is not a molblock
    mol_blocks.pop(0)

    if len(mol_blocks) > reactant_count + product_count:
        raise ValueError(
            "The number of mol blocks in the rxn block is greater than the number of reactants and products."  # noqa: E501
        )

    # if any element does not end with "M END", then remove all lines "M END".
    # This is important to capture the end of the last molblock correctly and
    # not include reaction data
    for ii, mol_block in enumerate(mol_blocks):
        if not mol_block.endswith(end_string):
            mol_blocks[ii] = mol_block.split(end_string)[0] + end_string  # add the end string back on

    # now we have the molblocks, we can split them into reactants and products
    reactants = mol_blocks[:reactant_count]
    products = mol_blocks[reactant_count:]

    return reactants, products
