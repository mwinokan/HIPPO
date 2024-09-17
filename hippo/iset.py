import mcol

import logging

logger = logging.getLogger("HIPPO")


class InteractionTable:
    """Class representing all :class:`.Interaction` objects in the 'interaction' table of the :class:`.Database`.

    .. attention::

            :class:`.InteractionTable` objects should not be created directly. Instead use the :meth:`.HIPPO.interactions` property.

    """

    def __init__(self, db: "Database", table: str = "interaction") -> None:

        self._db = db
        self._df = None
        self._table = table

    ### PROPERTIES

    @property
    def db(self) -> "Database":
        """Returns the associated :class:`.Database`"""
        return self._db

    @property
    def table(self) -> str:
        """Returns the name of the :class:`.Database` table"""
        return self._table

    @property
    def df(self) -> "pandas.DataFrame":
        """DataFrame representation of the interactions

        :returns: a ``pandas.Dataframe`` of the interactions

        """

        if self._df is None:
            records = self.db.select_all_where(
                table="interaction", key=f"interaction_id > 0", multiple=True
            )
            df = df_from_interaction_records(self.db, records)
            self._df = df

        return self._df

    ### DUNDERS

    def __len__(self) -> int:
        """The total number of interactions"""
        return self.db.count(self.table)


class InteractionSet:
    """Class representing a subset of the :class:`.Interaction` objects in the 'interaction' table of the :class:`.Database`.

    .. attention::

            :class:`.InteractionSet` objects should not be created directly. Instead use :meth:`.Pose.interactions`, or :meth:`.PoseSet.interactions` methods.

    """

    def __init__(
        self,
        db: "Database",
        indices: list = None,
        table: str = "interaction",
    ) -> None:

        self._db = db
        self._table = table

        indices = indices or []

        if not isinstance(indices, list):
            indices = list(indices)

        indices = [int(i) for i in indices]

        self._indices = sorted(list(set(indices)))
        self._df = None

    ### FACTORIES

    @classmethod
    def from_pose(
        cls, pose: "Pose | PoseSet", table: str = "interaction"
    ) -> "InteractionSet":
        """Construct a :class:`.InteractionSet` from one or more poses.

        :param pose: a :class:`.Pose` or :class:`.PoseSet` object
        :returns: an :class:`.InteractionSet`

        """

        self = cls.__new__(cls)

        ### get the ID's

        from .pset import PoseSet

        if isinstance(pose, PoseSet):

            # check if all poses have fingerprints
            (has_invalid_fps,) = pose.db.select_where(
                query="COUNT(1)",
                table="pose",
                key=f"pose_id IN {pose.str_ids} AND pose_fingerprint = 0",
            )

            if has_invalid_fps:
                logger.warning(f"{has_invalid_fps} Poses have not been fingerprinted")

            sql = f"""
            SELECT interaction_id FROM {table}
            WHERE interaction_pose IN {pose.str_ids}
            """

        else:

            sql = f"""
            SELECT interaction_id FROM {table}
            WHERE interaction_pose = {pose.id}
            """

        ids = pose.db.execute(sql).fetchall()

        ids = [i for i, in ids]

        self.__init__(pose.db, ids, table=table)

        return self

    @classmethod
    def from_residue(
        cls,
        db: "Database",
        residue_number: int,
        chain: None | str = None,
        target: "Target | int" = 1,
    ) -> "InteractionSet":
        """Get the set of interactions for a given residue number (and chain)

        :param db: HIPPO :class:`.Database`
        :param residue_number: the residue number
        :param chain: the chain name / letter, defaults to any chain
        :param target: the protein :class:`.Target` object or ID, defaults to first target in database
        :returns: a :class:`.InteractionSet` object
        """

        from .target import Target

        self = cls.__new__(cls)

        if isinstance(target, Target):
            target = target.id

        sql = f"""
        SELECT interaction_id FROM interaction
        INNER JOIN feature
        ON interaction_feature = feature_id
        WHERE feature_target = {target}
        AND feature_residue_number = {residue_number}
        """

        if chain:
            sql += f' AND feature_chain_name = "{chain}"'

        ids = db.execute(sql).fetchall()

        ids = [i for i, in ids]

        self.__init__(db, ids)

        return self

    ### PROPERTIES

    @property
    def indices(self) -> list[int]:
        """Returns the ids of interactions in this set"""
        return self._indices

    @property
    def ids(self) -> list[int]:
        """Returns the ids of interactions in this set"""
        return self._indices

    @property
    def db(self) -> "Database":
        """The associated HIPPO :class:`.Database`"""
        return self._db

    @property
    def table(self) -> str:
        """Get the name of the database table"""
        return self._table

    @property
    def str_ids(self) -> str:
        """Return an SQL formatted tuple string of the :class:`.Interaction` IDs"""
        return str(tuple(self.ids)).replace(",)", ")")

    @property
    def classic_fingerprint(self) -> dict:
        """Classic HIPPO fingerprint dictionary, mapping protein :class:`.Feature` ID's to the number of corresponding ligand features (from any :class:`.Pose`)"""
        return self.get_classic_fingerprint()

    @property
    def df(self) -> "pandas.DataFrame":
        """DataFrame representation of the interactions

        :returns: a ``pandas.Dataframe`` of the interactions

        """

        if self._df is None:
            records = self.db.select_all_where(
                table=self.table,
                key=f"interaction_id IN {self.str_ids}",
                multiple=True,
            )
            df = df_from_interaction_records(self.db, records)
            self._df = df

        return self._df

    @property
    def residue_number_chain_pairs(self) -> list[tuple[int]]:
        """Get a list of ``(residue_number, chain_name)`` tuples"""

        sql = f"""
        SELECT DISTINCT feature_residue_number, feature_chain_name FROM {self.table}
        INNER JOIN feature
        ON feature_id = interaction_feature
        WHERE interaction_id IN {self.str_ids}
        """

        return self.db.execute(sql).fetchall()

    @property
    def num_features(self) -> int:
        """Count the funmber of protein :class:`.Feature`s with which interactions are formed"""

        (count,) = self.db.execute(
            f"""
        SELECT COUNT(DISTINCT interaction_feature) FROM {self.table}
        WHERE interaction_id IN {self.str_ids}
        """
        ).fetchone()

        return count

    @property
    def avg_num_interactions_per_feature(self) -> float:
        """Average number of interactions formed with each protein :class:`.Feature`"""

        (count,) = self.db.execute(
            f"""
        WITH counts AS
        (
            SELECT interaction_feature, COUNT(1) AS count FROM {self.table}
            WHERE interaction_id IN {self.str_ids}
            GROUP BY interaction_feature
        )

        SELECT AVG(count) FROM counts
        """
        ).fetchone()

        return count

    @property
    def per_feature_count_std(self) -> float:
        """A measure for how evenly protein :class:`.Feature`s are being interacted with"""

        counts = self.db.execute(
            f"""
        SELECT interaction_feature, COUNT(1) AS count FROM interaction
        WHERE interaction_id IN {self.str_ids}
        GROUP BY interaction_feature
        """
        ).fetchone()

        counts = [c for c, in counts]

        from numpy import std

        return -std(counts)

    ### METHODS

    def summary(
        self,
        families: bool = False,
    ) -> None:
        """Print a summary of this :class:`.InteractionSet`"""

        logger.header(self)

        for interaction in self:
            # print(interaction)

            # logger.var(f'{interaction.family_str}', f'{interaction.distance:.1f}')
            s = f"{interaction.description}"

            if families:
                s += f" {interaction.feature.family} ~ {interaction.family}"

            logger.var(s, f"{interaction.distance:.1f}", dict(unit="Å"))

    def get_classic_fingerprint(self) -> dict:
        """Classic HIPPO fingerprint dictionary, mapping protein :class:`.Feature` ID's to the number of corresponding ligand features (from any :class:`.Pose`)"""

        pairs = self.db.execute(
            f"""
        SELECT interaction_feature, COUNT(1) FROM {self.table}
        WHERE interaction_id IN {self.str_ids}
        GROUP BY interaction_feature
        """
        ).fetchall()

        return {f: c for f, c in pairs}

    def resolve(
        self,
        debug: bool = False,
        commit: bool = True,
        # table: str = 'interaction',
    ) -> "InteractionSet":
        """Resolve into predicted key interactions. In place modification.

        :param debug: Increased verbosity for debugging (Default value = False)
        :param commit: commit the changes (Default value = True)
        :returns: a filtered :class:`.InteractionSet`
        """

        keep_list = []

        table = self.table

        ### H-Bonds (closest)

        sql = f"""
        SELECT interaction_id, MIN(interaction_distance)
        FROM {table}
        WHERE interaction_id IN {self.str_ids}
        AND interaction_type = "Hydrogen Bond"
        GROUP BY interaction_atom_ids
        """

        records = self.db.execute(sql).fetchall()
        ids = [a for a, b in records]
        keep_list += ids

        ### pi-stacking (closest)

        sql = f"""
        SELECT interaction_id, MIN(interaction_distance)
        FROM {table}
        INNER JOIN feature
        ON feature_id = interaction_feature
        WHERE interaction_id IN {self.str_ids}
        AND interaction_type = "π-stacking"
        GROUP BY feature_atom_names
        """
        # GROUP BY interaction_atom_ids

        records = self.db.execute(sql).fetchall()
        ids = [a for a, b in records]
        keep_list += ids

        ### pi-cation (closest)

        sql = f"""
        SELECT interaction_id, MIN(interaction_distance)
        FROM {table}
        WHERE interaction_id IN {self.str_ids}
        AND interaction_type = "π-cation"
        GROUP BY interaction_atom_ids
        """
        # GROUP BY interaction_atom_ids

        records = self.db.execute(sql).fetchall()
        ids = [a for a, b in records]
        keep_list += ids

        ### electrostatic (closest)

        sql = f"""
        SELECT interaction_id, MIN(interaction_distance)
        FROM {table}
        WHERE interaction_id IN {self.str_ids}
        AND interaction_type = "Electrostatic"
        GROUP BY interaction_atom_ids
        """
        # GROUP BY interaction_atom_ids

        records = self.db.execute(sql).fetchall()
        ids = [a for a, b in records]
        keep_list += ids

        ### hydrophobic

        sql = f"""
        SELECT interaction_id, interaction_distance
        FROM {table}
        WHERE interaction_id IN {self.str_ids}
        AND interaction_type = "Hydrophobic"
        """

        records = self.db.execute(sql).fetchall()
        ids = [a for a, b in records]
        subset = InteractionSet(self.db, ids, table=table)

        # aggregate lumped

        hydrophobic_interactions_in_lumped = {}
        lumped_hydrophobic_in_lumped_lumped = {}

        for interaction in subset:
            families = (interaction.feature.family, interaction.family)

            if families == ("LumpedHydrophobe", "Hydrophobe"):
                for name in interaction.feature.atom_names.split():
                    key = (name, interaction.atom_ids[0])
                    if key not in hydrophobic_interactions_in_lumped:
                        hydrophobic_interactions_in_lumped[key] = []
                    hydrophobic_interactions_in_lumped[key].append(interaction.id)

            elif families == ("Hydrophobe", "LumpedHydrophobe"):
                for atom_id in interaction.atom_ids:
                    key = (interaction.feature.atom_names, atom_id)
                    if key not in hydrophobic_interactions_in_lumped:
                        hydrophobic_interactions_in_lumped[key] = []
                    hydrophobic_interactions_in_lumped[key].append(interaction.id)

            elif families == ("LumpedHydrophobe", "LumpedHydrophobe"):
                for name in interaction.feature.atom_names.split():
                    for atom_id in interaction.atom_ids:
                        key = (name, atom_id)
                        if key not in hydrophobic_interactions_in_lumped:
                            hydrophobic_interactions_in_lumped[key] = []
                        hydrophobic_interactions_in_lumped[key].append(interaction.id)

                key = interaction.feature.atom_names
                lumped_hydrophobic_in_lumped_lumped[key] = tuple(interaction.atom_ids)

        keep_hydrophobic_ids = set(subset.ids)
        rev_hydrophobic_in_lumped_lumped = {
            v: k for k, v in lumped_hydrophobic_in_lumped_lumped.items()
        }

        # modify keep list by those covered in lumped

        for interaction in subset:
            families = (interaction.feature.family, interaction.family)

            if families == ("Hydrophobe", "Hydrophobe"):
                key = (interaction.feature.atom_names, interaction.atom_ids[0])

                if key in hydrophobic_interactions_in_lumped:
                    keep_hydrophobic_ids -= set([interaction.id])

            elif families == ("LumpedHydrophobe", "Hydrophobe"):

                key = interaction.feature.atom_names

                if key in lumped_hydrophobic_in_lumped_lumped:
                    atom_id = interaction.atom_ids[0]
                    value = lumped_hydrophobic_in_lumped_lumped[key]
                    if atom_id in value:
                        keep_hydrophobic_ids -= set([interaction.id])

            elif families == ("Hydrophobe", "LumpedHydrophobe"):

                key = tuple(interaction.atom_ids)

                if key in rev_hydrophobic_in_lumped_lumped:

                    atom_name = interaction.feature.atom_names
                    value = rev_hydrophobic_in_lumped_lumped[key]

                    if atom_name in value:
                        keep_hydrophobic_ids -= set([interaction.id])

        keep_list += list(keep_hydrophobic_ids)

        ### cull non-keepers

        cull_list = set(self.ids) - set(keep_list)
        cull_iset = InteractionSet(self.db, cull_list)
        self.db.delete_where(
            table=table,
            key=f"interaction_id IN {cull_iset.str_ids}",
            commit=commit,
        )
        self._indices = sorted(list(set(keep_list)))

        ### revisit hydrophobes

        # for a given protein feature, choose the closest interaction

        cull_list = []

        hydrophobic_keeper_iset = InteractionSet(self.db, keep_hydrophobic_ids)

        sql = f"""
        SELECT interaction_id, MIN(interaction_distance)
        FROM {table}
        WHERE interaction_id IN {hydrophobic_keeper_iset.str_ids}
        GROUP BY interaction_feature
        """

        records = self.db.execute(sql).fetchall()
        ids = [a for a, b in records]

        cull_list = set(hydrophobic_keeper_iset.ids) - set(ids)
        cull_iset = InteractionSet(self.db, cull_list)
        self.db.delete_where(
            table=table,
            key=f"interaction_id IN {cull_iset.str_ids}",
            commit=commit,
        )
        self._indices = sorted(list(set(keep_list) - cull_list))

        ### Summary

        if debug:
            self.summary()

    ### DUNDERS

    def __len__(self) -> int:
        """The number of interactions in this set"""
        return len(self.indices)

    def __repr__(self) -> str:
        """Formatted command-line representation"""
        return (
            f"{mcol.bold}{mcol.underline}"
            "{"
            f"I x {len(self)}"
            "}"
            f"{mcol.unbold}{mcol.ununderline}"
        )

    def __iter__(self):
        """Iterate through interactions in this set"""
        return iter(
            self.db.get_interaction(id=i, table=self.table) for i in self.indices
        )

    def __getitem__(self, key) -> "Interaction | InteractionSet":
        """Get interaction or subsets thereof from this set"""
        match key:
            case int():
                index = self.indices[key]
                return self.db.get_interaction(id=index, table=self.table)

            case slice():
                indices = self.indices[key]
                return InteractionSet(self.db, indices, table=self.table)

            case _:
                raise NotImplementedError


def df_from_interaction_records(
    db: "Database",
    records: list[tuple],
) -> "pandas.DataFrame":
    """Construct a dataframe from the 'interaction' table records"""

    import json
    from pandas import DataFrame

    data = []
    for record in records:

        (
            id,
            feature_id,
            pose_id,
            type,
            family,
            atom_ids,
            prot_coord,
            lig_coord,
            distance,
            angle,
            energy,
        ) = record

        feature = db.get_feature(id=feature_id)

        d = dict(id=id)

        d["feature_id"] = feature_id
        d["pose_id"] = pose_id
        d["target_id"] = feature.target

        # d['type'] = INTERACTION_TYPES[(feature.family, family)]
        d["type"] = type

        d["prot_family"] = feature.family
        d["lig_family"] = family

        d["residue_name"] = feature.residue_name
        d["residue_number"] = feature.residue_number
        d["chain_name"] = feature.chain_name

        d["distance"] = distance
        d["angle"] = angle
        d["energy"] = energy

        d["prot_coord"] = json.loads(prot_coord)
        d["lig_coord"] = json.loads(lig_coord)

        d["prot_atoms"] = feature.atom_names
        d["lig_atoms"] = atom_ids

        data.append(d)

    df = DataFrame.from_records(data=data)

    return df
