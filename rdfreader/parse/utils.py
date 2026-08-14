import datetime
import logging
import re
from typing import Any, Callable, Optional

logger = logging.getLogger(__name__)

# http://c4.cabrillo.edu/404/ctfile.pdf
CTF_DEFAULT_MOLBLOCK_HEADER_FORMAT_STRING: str = "IIPPPPPPPPMMDDYYHHmmddSSssssssssssEEEEEEEEEEEERRRRRR"
CTF_RXNBLOCK_HEADER_FORMAT_STRING: str = "IIIIIIPPPPPPPPPMMDDYYYYHHmmRRRRRRR"
SPRESI_RXNBLOCK_HEADER_FORMAT_STRING: str = "IIIIIIPPPPPPPPPPMMDDYYHHmmRRRRRRR"
CTF_COMPONENT_COUNT_FORMAT_STRING: str = "rrrppp"
# default mapping of letter to field. This is used to map the letter in the
# mol block header to the field in the metadata.
# key: letter, value: tuple containing the field name, the field type and a
# default value (if the field is empty). Field typle should be a castable data
# type.
# If the field is empty, the default value will be used.
CTF_DEFAULT_LETTER_TO_FIELD_MAPPING: dict[str, tuple[str, Callable, Any]] = {
    "I": ("user_initials", str, ""),
    "P": ("program_name", str, ""),
    "M": ("month", int, 0),
    "D": ("day", int, 0),
    "Y": ("year", int, 0),
    "H": ("hour", int, 0),
    "m": ("minute", int, 0),
    "d": ("dimensional_codes", str, ""),
    "S": ("scaling_factor_1", int, 0),
    "s": ("scaling_factor_2", float, 0.0),
    "E": ("energy", float, 0.0),
    "R": ("registry_number", str, ""),
    "r": ("reactant_count", int, 0),
    "p": ("product_count", int, 0),
}


def _default_line_item(cast_type: Callable = str, default: Any = None) -> Any:
    """Return a default value for a line item."""
    if default is not None:
        return default

    if cast_type in (int, float):
        return cast_type(0)
    if cast_type == str:
        return ""
    else:
        return default


def get_line_item(
    line: str,
    character_index: Optional[tuple[Optional[int], Optional[int]]] = None,
    cast_type: Callable = str,
    default: Any = None,
    catch_casting_exceptions: bool = False,
) -> Any:
    """Process a line from an rdf/mol file or block.

    If character_index is provided, then the characters within that slice will
    be returned.
    If cast_type is provided, then the characters will be cast to that type.
    If the line is empty a default will be provided. For numeric types (as set
    by cast_type), the default will be 0, for strings it will be an empty
    string. For other types, the default will be None. This can be overridden
    by providing a default.
    If catch_casting_exceptions is True, then any casting exceptions will be
    caught and the default will be returned.

    Strips whitespace and newlines.

    Parameters
    ----------
    line : str
        A line.
    character_index : tuple[int, int]
        A tuple of the start and end character index.

    Returns
    -------
    str
        The line item.
    """
    line = line[slice(*character_index)] if character_index else line
    line = line.strip("\n")  # remove newline
    line = line.strip()  # remove whitespace

    if not line:
        # return a default value if the line is empty
        return _default_line_item(cast_type, default)

    if cast_type:
        try:
            # attempt to cast the line to the specified type
            line = cast_type(line)
        except ValueError:
            if catch_casting_exceptions:
                # return a default value if the line cannot be cast
                return _default_line_item(cast_type, default)
            else:
                raise

    return line


def get_whole_line_item(line: str) -> str:
    """Return a data item that is from a whole line.

    Item must be a string.
        Parameters
        ----------
        line : str
            A line.

        Returns
        -------
        str
            The data item.
    """
    return get_line_item(line, (0, len(line)), str, "")


def parse_format_string(format_string) -> dict[str, tuple]:
    """Parses a CTF format string (http://c4.cabrillo.edu/404/ctfile.pdf) and
    returns a dictionary where the key is the letter in the format string and
    the value is a tuple of the start and end character index.

    Parameters
    ----------
    format_string : str
        A CTF format string.

    Returns
    -------
    dict[str, tuple]
        A dictionary of the letter in the format string and the start and end
        character index as a tuple.
    """
    format_string_dict = {}

    last_letter = format_string[0]
    current_letter_start_index = 0  # the index where the current letter started
    for ii, letter in enumerate(format_string):
        if letter != last_letter:
            # if the letter has changed add the previous letter to the
            # dictionary
            format_string_dict[last_letter] = (current_letter_start_index, ii)
            current_letter_start_index = ii

        last_letter = letter

    # add the last letter to the dictionary
    format_string_dict[last_letter] = (current_letter_start_index, ii + 1)

    return format_string_dict


def dict_elements_to_datetime(
    dd: dict[str, Any],
    date_time_key: str = "date_time",
    delete_initial_keys: bool = True,
    catch_datetime_exceptions: bool = True,
) -> dict[str, Any]:
    """Search a dictionary for datetime elements and adds a a new key to the
    dictionary with a datetime object.

    Searches for datetime keys called: "hour", "minute", "second", "day",
    "month", "year" and adds a new key called <date_time_key> with a datetime
    object.

    If delete_initial_keys is set, the original keys will be deleted.

    If the datetime is not parseable, None will be added to the dictionary.

    Parameters
    ----------
    dd : dict[str, Any]
        A dictionary.
    date_time_key : str
        The key to add to the dictionary.
    delete_initial_keys : bool
        If True, the keys "hour", "minute", "second", "day", "month", "year"
        will be deleted from the dictionary.
    catch_datetime_exceptions
        If True, parsing exceptions from datetime.datetime will be caught and
        None will be added to the dictionary.

    Returns
    -------
    dict[str, Any]
        A dictionary with the new key added.
    """
    date_time_args = {k: dd[k] for k in ["hour", "minute", "second", "day", "month", "year"] if k in dd}

    if not date_time_args:
        # if no datetime keys are found, return the dictionary as is
        return dd

    try:
        dd[date_time_key] = datetime.datetime(**date_time_args)
    except (ValueError, TypeError):
        if catch_datetime_exceptions:
            logger.warning(f"Could not parse datetime from {dd}")
            dd[date_time_key] = None
        else:
            raise

    if delete_initial_keys:
        [dd.pop(key, None) for key in ["hour", "minute", "second", "day", "month", "year"]]

    return dd


def _parse_block_header_line(
    header_line: str,
    header_format_string: str,
    header_field_mapping: dict[str, tuple[str, Callable, Any]],
) -> dict[str, Any]:
    """Parse the header line of a rxn or mol block.

    Parameters
    ----------
    header_line : list[str]
        A single line from the header.
    header_format_string : str
        See get_mol_block_metadata for more information.
    header_field_mapping : dict[str, tuple(str, Callable, Any)]
        See get_mol_block_metadata for more information.


    Returns
    -------
    dict
        A dictionary of metadata.
    """
    format_string_dict = parse_format_string(header_format_string)

    metadata = {}
    for letter, character_index in format_string_dict.items():
        if letter in header_field_mapping:
            field_name = header_field_mapping[letter][0]
            data_type = header_field_mapping[letter][1]
            try:
                default_value = header_field_mapping[letter][2]
            except IndexError:
                default_value = None
            metadata[field_name] = get_line_item(header_line, character_index, data_type, default_value)
        else:
            raise ValueError(f"The letter {letter} does not appear in the format field mapping.")  # noqa: E501

    metadata = dict_elements_to_datetime(metadata)

    return metadata


def make_string_python_safe(string: str) -> str:
    """Remove/replace characters that are not allowed in python
    variable/function names.

    Parameters
    ----------
    string : str
        A string.

    Returns
    -------
    str
        The python safe string.
    """
    if string is None:
        return None

    string = string.strip()
    string = string.strip("\n")
    string = string.strip("\r")

    if string == "":
        return None

    # remove existing underscores
    # string = string.replace("_", "")

    # add a leading underscore if string starts with a number
    if string[0].isdigit():
        string = f"_{string}"

    # replace all non-alphanumeric characters with an underscore
    string = re.sub(r"[^a-zA-Z0-9_]", "_", string)

    # replace multiple underscores with a single underscore
    string = re.sub(r"_{2,}", "_", string)

    # lower case the string
    string = string.lower()

    # remove trailing underscores
    string = string.rstrip("_")

    return string


def parse_yield(yield_string: str) -> float:
    """Attempts to parse the yield string into a float.

    If muliple numbers are detected, will return the average.

    If no numbers are detected, will return None.

    Parameters
    ----------
    yield_string : str

    Returns
    -------
    float
        Yield as a float.
    """
    re_patterns = [
        r"^([0-9]+\.?[0-9]?)$",  # matches a single int or float
        r"^([0-9]+\.?[0-9]?)\s{0,}[-,;:]{0,}\s{0,}([0-9]+\.?[0-9]?)$",  # matches two ints or floats separated by a dash, comma, semicolon, colon or space  # noqa: E501
    ]

    for re_pattern in re_patterns:
        match = re.search(re_pattern, yield_string)
        if match:
            # get the match groups as a list
            # (first match group is the whole string)
            match_groups = match.groups()
            # convert the match groups to floats
            match_groups = [float(group) for group in match_groups]
            # average the match groups
            return sum(match_groups) / len(match_groups)

    # if we get here, then we didn't find a match
    logger.warning(f"Could not parse yield from '{yield_string}'. Returning None.")

    return None
