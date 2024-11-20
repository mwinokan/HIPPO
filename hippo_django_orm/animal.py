import mrich
from .database import Database


class HIPPO:
    """
    HIPPO class
    """

    def __init__(self, name: str, db_path: str) -> None:

        self._name = name
        self._db_path = db_path

        # DB must be initialised before importing any models
        self._db = Database(path=db_path)

        from .compound import CompoundTable
        from .pose import PoseTable
        from .target import TargetTable

        self.compounds = CompoundTable()
        self.poses = PoseTable()
        self.targets = TargetTable()

    ### FACTORIES

    ### PROPERTIES

    @property
    def name(self):
        return self._name

    @property
    def name(self) -> str:
        """Returns the project name

        :returns: project name
        """
        return self._name

    @property
    def db_path(self) -> str:
        """Returns the database path"""
        return self._db_path

    @property
    def db(self) -> Database:
        """Returns the Database object"""
        return self._db

    @property
    def num_compounds(self) -> int:
        """Total number of Compounds in the Database"""
        return len(self.compounds)

    @property
    def num_poses(self) -> int:
        """Total number of poses in the Database"""
        return len(self.poses)

    @property
    def num_targets(self) -> int:
        """Total number of targets in the Database"""
        return len(self.targets)

    ### METHODS

    def summary(self) -> None:
        """Print a text summary of this HIPPO"""
        mrich.header(self)
        mrich.var("db_path", self.db_path)
        mrich.var("#compounds", self.num_compounds)
        mrich.var("#poses", self.num_poses)
        mrich.var("#targets", self.num_targets)
        # mrich.var("#reactions", self.num_reactions)
        # mrich.var("#tags", self.num_tags)
        # mrich.var("tags", self.tags.unique)

    ### DUNDERS

    def __str__(self) -> str:
        """Unformatted string representation of this HIPPO"""
        return f'HIPPO("{self.name}")'

    def __repr__(self) -> str:
        """Returns a command line representation of this HIPPO"""
        return f"{mcol.bold}{mcol.underline}{self}{mcol.clear}"

    def __rich__(self) -> str:
        """Representation for mrich"""
        return f"[bold underline]{self}"
