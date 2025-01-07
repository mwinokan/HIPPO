import mcol
import mrich

from dataclasses import dataclass

from .target import Target


@dataclass
class Feature:
    """Pharmocophoric feature in a protein

    attributes:
            id: Database ID
            family: Feature family
            chain_name: Protein chain name/letter
            residue_name: Protein residue name
            residue_number: Protein residue number
            atom_names: Protein atom names (whitespace-delimited)

    """

    id: int
    family: Target
    target: int
    chain_name: str
    residue_name: str
    residue_number: int
    atom_names: str

    def __str__(self) -> str:
        """Unformatted string representation"""
        return f"F{self.id}"

    def __repr__(self) -> str:
        """ANSI Formatted string representation"""
        return f"{mcol.bold}{mcol.underline}{self.family} {self.chain_name} {self.residue_name} {self.residue_number} [{self.atom_names}]{mcol.unbold}{mcol.ununderline}"

    def __rich__(self) -> str:
        """Representation for mrich"""
        return f"[bold underline]{self.family} {self.chain_name} {self.residue_name} {self.residue_number} [{self.atom_names}]"

    @property
    def chain_res_name_number_str(self) -> str:
        """Return a string representation of the feature"""
        return f"{self.chain_name} {self.residue_name} {self.residue_number}"

    @property
    def res_name_number_str(self) -> str:
        """Return a string representation of the feature"""
        return f"{self.residue_name} {self.residue_number}"

    @property
    def res_number_name_tuple(self) -> str:
        """Return a tuple representation of the feature"""
        return (self.residue_number, self.residue_name)

    @property
    def res_name_number_family_str(self) -> str:
        return f"{self.residue_name} {self.residue_number} {self.family}"

    @property
    def backbone(self) -> bool:
        """Are any of the atoms referenced by this feature on the backbone?"""
        from molparse.amino import BB_NAMES

        for atom_name in self.atom_names.split(","):
            if atom_name in BB_NAMES:
                return True

        return False

    @property
    def sidechain(self) -> bool:
        """Are any of the atoms referenced by this feature on the sidechain?"""
        from molparse.amino import BB_NAMES

        for atom_name in self.atom_names.split(","):
            if atom_name not in BB_NAMES:
                return True

        return False
