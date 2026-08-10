import logging
from typing import Any, Callable

from rdfreader.parse.utils import (
    CTF_DEFAULT_LETTER_TO_FIELD_MAPPING,
    CTF_DEFAULT_MOLBLOCK_HEADER_FORMAT_STRING,
    _parse_block_header_line,
    get_line_item,
    get_whole_line_item,
)

logger = logging.getLogger(__name__)


def get_mol_block_metadata(
    mol_block: str,
    header_format_string: str = CTF_DEFAULT_MOLBLOCK_HEADER_FORMAT_STRING,
    header_field_mapping: dict[str, tuple[str, Callable, Any]] = CTF_DEFAULT_LETTER_TO_FIELD_MAPPING,
) -> dict[str, Any]:
    """Extract metadata from a mol block string.

    Parameters
    ----------
    mol_block : str
        A mol block string.
    header_format_string : str
        The format string for the line 2 of the mol_block header. Default
        value is according to the mol_block spec at
        http://c4.cabrillo.edu/404/ctfile.pdf. The letters must appear in the
        spec, but the order can be changed.
    header_field_mapping : dict[str, tuple(str, Callable, Any)]
        A dictionary which maps a letter in the header format string to a
        field name. The field name is the key and the letter is the value. The
        default mapping is according to the mol_block spec at
        http://c4.cabrillo.edu/404/ctfile.pdf.
        key: letter, value: tuple containing the field name, the field type
        and a default value (if the field is empty). Field typle should be a
        castable data type.
        If the field is empty, the default value will be used.

    Returns
    -------
    dict
        A dictionary of metadata.
    """
    metadata = {}
    mol_block_lines = mol_block.split("\n")
    metadata["molecule_name"] = get_whole_line_item(mol_block_lines[0])
    metadata.update(_parse_block_header_line(mol_block_lines[1], header_format_string, header_field_mapping))
    metadata["comment"] = get_whole_line_item(mol_block_lines[2])
    large_regno: str = _parse_large_regno(mol_block_lines)
    if large_regno is not None:
        # overwrite the registry number from the header with the large
        # registry number if it is present
        metadata["registry_number"] = large_regno

    return metadata


def _parse_large_regno(mol_block_lines: list[str]) -> str:
    """Searches the molblock a line beginning with M REG and returns the value
    of it if present.

    Parameters
    ----------
    mol_block_lines : list[str]
        A molblock string that has been split into a list of lines.

    Returns
    -------
    Optional[str]
        The registry number. If None is returned, the registry number is not
        present.
    """
    for line in mol_block_lines:
        if line.startswith("M  REG "):
            return get_line_item(line, [7, len(line)])

    return None
