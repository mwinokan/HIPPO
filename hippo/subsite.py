import mcol


class Subsite:
    """
    Class representing a subsite/subsite on a protein :class:`.Target`

    .. attention::

        :class:`.Subsite` objects should not be created directly. Instead use :meth:`.Target.subsites`.

    """

    _table = "subsite"

    def __init__(self, db: "Database", id: int, target_id: int, name: str):

        self._db = db
        self._id = id
        self._target_id = target_id
        self._name = name
        self._metadata = None
        self._target = None

    ### FACTORIES

    ### PROPERTIES

    @property
    def db(self) -> "Database":
        """Returns a pointer to the parent database"""
        return self._db

    @property
    def id(self) -> int:
        """Returns the SubsiteTag's database ID"""
        return self._id

    @property
    def table(self):
        """Returns the name of the :class:`.Database` table"""
        return self._table

    @property
    def target(self) -> "Target":
        """Returns the associated protein :class:`.Target`"""
        if self._target is None:
            self._target = self.db.get_target(id=self.target_id)
        return self._target

    @property
    def target_id(self) -> int:
        """Returns the associated protein :class:`.Target` ID"""
        return self._target_id

    @property
    def name(self) -> str:
        """The subsite name"""
        return self._name

    @property
    def metadata(self) -> "MetaData":
        """Returns the SubsiteTag's metadata"""
        if self._metadata is None:
            self._metadata = self.db.get_metadata(table="subsite", id=self.id)
        return self._metadata

    @property
    def poses(self) -> "PoseSet | None":
        """Return all poses in this subsite"""
        from .pset import PoseSet

        indices = self.db.select_where(
            table="subsite_tag",
            query="subsite_tag_pose",
            multiple=True,
            key="ref",
            value=self.id,
        )
        indices = [i for i, in indices]
        if not indices:
            return None
        return PoseSet(self.db, indices, name="poses in {self}")

    ### METHODS

    ### DUNDERS

    def __str__(self):
        """Unformatted string representation"""
        return f"{self.target.name}->{self.name}"

    def __repr__(self):
        """Formatted string representation"""
        return f'{mcol.bold}{mcol.underline}S{self.id} "{self}"{mcol.unbold}{mcol.ununderline}'


class SubsiteTag:
    """
    Class representing a tag assigning a :class:`.Subsite` to a :class:`.Pose`

    .. attention::

        :class:`.SubsiteTag` objects should not be created directly. Instead use :meth:`.Pose.subsites`.

    """

    _table = "subsite_tag"

    def __init__(self, db: "Database", id: int, subsite_id: int, pose_id: int):

        self._db = db
        self._id = id
        self._subsite_id = subsite_id
        self._pose_id = pose_id
        self._metadata = None

        self._name = db.get_subsite_name(id=subsite_id)

    ### FACTORIES

    ### PROPERTIES

    @property
    def db(self) -> "Database":
        """Returns a pointer to the parent database"""
        return self._db

    @property
    def id(self) -> int:
        """Returns the SubsiteTag's database ID"""
        return self._id

    @property
    def table(self):
        """Returns the name of the :class:`.Database` table"""
        return self._table

    @property
    def pose_id(self) -> int:
        """Returns the associated :class:`.Pose`'s database ID"""
        return self._pose_id

    @property
    def subsite_id(self):
        """Returns the associated :class:`.Subsite`'s database ID"""
        return self._subsite_id

    @property
    def name(self):
        """Returns the associated :class:`.Subsite`'s name"""
        return self._name

    @property
    def metadata(self) -> "MetaData":
        """Returns the SubsiteTag's metadata"""
        if self._metadata is None:
            self._metadata = self.db.get_metadata(table="subsite_tag", id=self.id)
        return self._metadata

    ### METHODS

    ### DUNDERS

    def __str__(self):
        """Unformatted string representation"""
        return self.name

    def __repr__(self):
        """Formatted string representation"""
        return f'{mcol.bold}{mcol.underline}"{self}"{mcol.unbold}{mcol.ununderline}'
