import mcol

import logging

logger = logging.getLogger("HIPPO")


class Interaction:
    """A :class:`.Interaction` represents an interaction between an rdkit Feature on a :class:`.Pose` and a :class:`.Feature` on the protein :class:`.Target`.

    .. attention::

            :class:`.Interaction` objects should not be created directly. Instead use :meth:`.Pose.interactions`, or :meth:`.PoseSet.interactions` methods.

    """

    def __init__(
        self,
        db: "Database",
        id: int,
        feature_id: int,
        pose_id: int,
        type: str,
        family: str,
        atom_ids: str,
        prot_coord: str,
        lig_coord: str,
        distance: float,
        angle: float,
        energy: float | None,
        table: str = "interaction",
    ) -> None:

        import json

        # from interaction table
        self._id = id
        self._feature_id = feature_id
        self._pose_id = pose_id
        self._type = type
        self._family = family
        self._atom_ids = json.loads(atom_ids)
        self._prot_coord = json.loads(prot_coord)
        self._lig_coord = json.loads(lig_coord)
        self._distance = distance
        self._angle = angle
        self._energy = energy

        # placeholders
        self._pose = None
        self._feature = None
        self._table = table

        self._db = db

    ### PROPERTIES

    @property
    def id(self) -> int:
        """Returns the interaction's database ID"""
        return self._id

    @property
    def table(self) -> str:
        """Returns the name of the :class:`.Database` table"""
        return self._table

    @property
    def db(self) -> "Database":
        """Returns a pointer to the parent database"""
        return self._db

    @property
    def family(self) -> str:
        """The Feature family"""
        return self._family

    @property
    def pose_id(self) -> int:
        """Returns the associated :class:`.Pose`'s database ID"""
        return self._pose_id

    @property
    def pose(self) -> "Pose":
        """Returns the associated :class:`.Pose`'s object"""
        if not self._pose:
            self._pose = self.db.get_pose(id=self.pose_id)
        return self._pose

    @property
    def feature_id(self) -> int:
        """Returns the associated :class:`.Feature`'s database ID"""
        return self._feature_id

    @property
    def feature(self) -> "Feature":
        """Returns the associated :class:`.Feature`'s object"""
        if not self._feature:
            self._feature = self.db.get_feature(id=self.feature_id)
        return self._feature

    @property
    def atom_ids(self) -> list[int]:
        """Returns the indices of atoms making up the ligand feature"""
        return self._atom_ids

    @property
    def prot_coord(self) -> list[float]:
        """Returns the cartesian position of the protein :class:`.Feature`"""
        return self._prot_coord

    @property
    def lig_coord(self) -> list[float]:
        """Returns the cartesian position of the ligand feature"""
        return self._lig_coord

    @property
    def distance(self) -> float:
        """Returns the euclidian distance of the interaction"""
        return self._distance

    @property
    def angle(self) -> float | None:
        """Returns the interaction angle (only defined for π-stacking and π-cation interactions)"""
        return self._angle

    @property
    def energy(self) -> float | None:
        """Returns the interaction energy, if defined"""
        return self._energy

    @property
    def family_str(self) -> str:
        """String of the two feature families"""
        return f"{repr(self.feature)} ~ {self.family}"

    @property
    def type(self) -> str:
        """Interaction type string"""
        from molparse.rdkit.features import INTERACTION_TYPES

        return INTERACTION_TYPES[(self.feature.family, self.family)]

    @property
    def description(self) -> str:
        """One line description of this interaction"""
        s = f"{self.type} [{self.feature.chain_res_name_number_str}] {self.distance:.1f} Å"
        if self.angle:
            s += f", {self.angle:.1f} degrees"
        return s

    ### METHODS

    def summary(self) -> None:
        """Print a summary of this interaction's properties"""

        logger.header(f"Interaction {self.id}")

        logger.var("feature", self.feature)
        logger.var("pose", self.pose)
        logger.var("family", self.family)
        logger.var("atom_ids", self.atom_ids)
        logger.var("prot_coord", self.prot_coord)
        logger.var("lig_coord", self.lig_coord)
        logger.var("distance", self.distance)
        logger.var("angle", self.angle)
        logger.var("energy", self.energy)

    ### DUNDERS

    def __str__(self) -> str:
        """Plain string representation"""
        return f"I{self.id}"

    def __repr__(self) -> str:
        """Formatted string representation"""
        return f"{mcol.bold}{mcol.underline}{str(self)}{mcol.unbold}{mcol.ununderline}"
