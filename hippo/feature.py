import mcol

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
        """Formatted string representation"""
        return f"{mcol.bold}{mcol.underline}{self.family} {self.chain_name} {self.residue_name} {self.residue_number} [{self.atom_names}]{mcol.unbold}{mcol.ununderline}"

    @property
    def chain_res_name_number_str(self) -> str:
        """Return a string representation of the feature"""
        return f"{self.chain_name} {self.residue_name} {self.residue_number}"

    @property
    def res_name_number_str(self) -> str:
        """Return a string representation of the feature"""
        return f"{self.residue_name} {self.residue_number}"

    @property
    def res_name_number_family_str(self) -> str:
        return f"{self.residue_name} {self.residue_number} {self.family}"
