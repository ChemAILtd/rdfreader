# typing
from io import TextIOWrapper
from pathlib import Path
from typing import Optional

from rdfreader.chem.reaction import Reaction
from rdfreader.parse.utils import CTF_RXNBLOCK_HEADER_FORMAT_STRING


def parse_rdf_reg_num(line: str):
    return line.replace("$RFMT $RIREG ", "").strip()


class RDFParser:
    _header_retrieved: bool = False
    lineno: int = 1
    rdf_metadata: dict[str, str] = {}
    _next_line: Optional[str] = None  # Lookahead buffer

    def __init__(
        self,
        f: TextIOWrapper,
        header_format_string: str = CTF_RXNBLOCK_HEADER_FORMAT_STRING,
        except_on_invalid_molecule: bool = True,
        except_on_invalid_reaction: bool = True,
        parse_conditions: bool = True,
    ):
        """
        Parameters
        ----------
        f : TextIOWrapper
            The file to parse.
        header_format_string : str, optional
            The format string to use to parse the header of each rxn block.
        except_on_invalid_molecule : bool, optional
            If True, raise an exception if a molecule is invalid.
        except_on_invalid_reaction : bool, optional
            If True, raise an exception if a reaction is invalid.
        rdf_file_name : str, optional
            The name of the rdf file.
        parse_conditions : bool, optional
            Whether to parse the conditions of the reaction.
        """

        self.f = f
        self.header_format_string = header_format_string
        self.except_on_invalid_molecule = except_on_invalid_molecule
        self.except_on_invalid_reaction = except_on_invalid_reaction
        self.rdf_file_name = Path(f.name).name
        self.parse_conditions = parse_conditions

    def __iter__(self):
        return self

    def __next__(self):
        return self.next_reaction()

    def next_reaction(self):
        """Returns the next rxn block from a rdf file.

        Parameters
        ----------
        f : TextIOWrapper
            The file to parse.

        Returns
        -------
        tuple[str, str, int]
            The rxn block, the rxn id, and the line number of the start of the
            rxn block.
        """

        # get the next rxn block
        reaction = None
        rxn_block, rxn_id, start_lineno = self._next_rxn_block()
        try:
            reaction = Reaction(
                rxn_block=rxn_block,
                id=rxn_id,
                rdf_metadata=self.rdf_metadata,
                header_format_string=self.header_format_string,
                except_on_invalid_molecule=self.except_on_invalid_molecule,
                lineno=start_lineno,
                rdf_file=self.rdf_file_name,
            )

        except Exception as e:
            if self.except_on_invalid_reaction:
                raise e

        return reaction

    def _header(self):
        """Parse the header of a RDF file.

        Parameters
        ----------
        f : TextIOWrapper
            The file to parse.

        Returns
        -------
        dict[str, str]
            The version and date of the RDF file. These are just treated as
            strings as they are typically ignored and the structure of the
            datetime field is not defined in the specification.
        """
        if not self._header_retrieved:
            self._header_retrieved = True
            version = self._readline()[8:].strip()
            date_stamp = self._readline()[6:].strip()
            self.rdf_metadata = {"version": version, "date_stamp": date_stamp}

    def _next_rxn_block(self) -> tuple[str, str, int]:
        """Returns the next rxn block from a rdf file. If the end.

        Returns
        -------
        tuple[str, str, int]
            The rxn block, the rxn id, and the line number of the start of the
            rxn block.
        """

        self._header()
        start_lineno: int = self.lineno
        line: str = self._readline()
        if line == "":
            raise StopIteration
        # parse the rxn block deliminators
        if line.startswith("$RFMT"):
            # capture the reg number
            reg_no: str = parse_rdf_reg_num(line)
            line: str = self._readline()
        else:
            #  there is a problem with the file format, raise an exception
            raise Exception(f"Invalid RDF file format. Expected $RFMT, got {line} " f"at line {self.lineno}")

        # capture the rxn block
        rxn_block: str = ""
        while not line.startswith("$RFMT") and not line == "":
            rxn_block += line
            line = self._readline()

        # Store lookahead line for next iteration
        self._next_line = line
        self.lineno -= 1

        return rxn_block, reg_no, start_lineno

    def _readline(self):
        """Wraps f.read and increments the line number."""
        self.lineno += 1

        # Check the lookahead buffer first
        if self._next_line is not None:
            line = self._next_line
            self._next_line = None
            return line

        return self.f.readline()
