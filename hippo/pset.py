from .db import Database
from .pose import Pose
from .cset import IngredientSet

import mcol

from typing import Callable

import os

import mrich


class PoseTable:
    """Class representing all :class:`.Pose` objects in the 'pose' table of the :class:`.Database`.

    .. attention::

            :class:`.PoseTable` objects should not be created directly. Instead use the :meth:`.HIPPO.poses` property. See :doc:`getting_started` and :doc:`insert_elaborations`.

    Use as an iterable
    ==================

    Iterate through :class:`.Pose` objects in the table:

    ::

            for pose in animal.poses:
                ...


    Selecting poses in the table
    ============================

    The :class:`.PoseTable` can be indexed with :class:`.Pose` IDs, names, aliases, or list/sets/tuples/slices thereof:

    ::

            ptable = animal.poses

            # indexing individual compounds
            pose = ptable[13]                            # using the ID
            pose = ptable["BSYNRYMUTXBXSQ-UHFFFAOYSA-N"] # using the InChIKey
            pose = ptable["Ax0310a"]                     # using the alias

            # getting a subset of compounds
            pset = ptable[13,15,18]      # using IDs (tuple)
            pset = ptable[[13,15,18]]    # using IDs (list)
            pset = ptable[set(13,15,18)] # using IDs (set)
            pset = ptable[13:18]         # using a slice

    Tags and target IDs can also be used to filter:

    ::

            pset = animal.poses(tag='hits') # select compounds tagged with 'hits'
            pset = animal.poses(target=1)   # select poses from the first target

    """

    _table = "pose"
    _name = "all poses"

    def __init__(
        self,
        db: Database,
    ) -> None:

        self._db = db

    ### FACTORIES

    ### PROPERTIES

    @property
    def db(self) -> Database:
        """Returns the associated :class:`.Database`"""
        return self._db

    @property
    def table(self) -> str:
        """Returns the name of the :class:`.Database` table"""
        return self._table

    @property
    def name(self) -> str | None:
        """Returns the name of set"""
        return self._name

    @property
    def names(self) -> list[str]:
        """Returns the aliases of child poses"""
        return [p.name for p in self]

    @property
    def aliases(self) -> list[str]:
        """Returns the aliases of child poses"""
        result = self.db.select(table=self.table, query="pose_alias", multiple=True)
        return [q for q, in result]

    @property
    def inchikeys(self) -> list[str]:
        """Returns the inchikeys of child poses"""
        result = self.db.select(table=self.table, query="pose_inchikey", multiple=True)
        return [q for q, in result]

    @property
    def ids(self) -> list[int]:
        """Returns the IDs of child poses"""
        result = self.db.select(table=self.table, query="pose_id", multiple=True)
        return [q for q, in result]

    @property
    def tags(self) -> set[str]:
        """Returns the set of unique tags present in this pose set"""
        values = self.db.select_where(
            table="tag",
            query="DISTINCT tag_name",
            key="tag_pose IS NOT NULL",
            multiple=True,
        )
        return set(v for v, in values)

    @property
    def num_fingerprinted(self) -> int:
        """Count the number of fingerprinted poses"""
        return self.db.count_where(
            table="pose",
            key="fingerprint",
            value=1,
        )

    @property
    def id_name_dict(self) -> dict[int, str]:
        """Return a dictionary mapping pose ID's to their name"""

        records = self.db.select(
            table=self.table, query="pose_id, pose_inchikey, pose_alias", multiple=True
        )

        lookup = {}
        for i, inchikey, alias in records:
            if alias:
                lookup[i] = alias
            else:
                lookup[i] = inchikey

        return lookup

    ### METHODS

    def get_by_tag(
        self,
        tag: str,
        inverse: bool = False,
    ) -> "PoseSet":
        """Get all child poses with a certain tag

        :param tag: tag to search for
        :param inverse: invert the selection
        :returns: a :class:`.PoseSet` of the subset

        """

        if not inverse:

            values = self.db.select_where(
                query="tag_pose", table="tag", key="name", value=tag, multiple=True
            )

        else:

            values = self.db.select_where(
                query="tag_pose", table="tag", key="name", value=tag, multiple=True
            )

            if not values:
                return self

            ids = [v for v, in values if v]

            values = self.db.select_where(
                query="pose_id",
                table="pose",
                key=f"pose_id NOT IN {self.str_ids}",
                multiple=True,
            )

        if not values:
            return None

        ids = [v for v, in values if v]

        pset = self[ids]

        if inverse:
            pset._name = f'poses not tagged "{tag}"'
        else:
            pset._name = f'poses tagged "{tag}"'
        return pset

    def get_by_target(
        self,
        *,
        id: int,
    ) -> "PoseSet":
        """Get all child poses with a certain :class:`.Target` ID:

        :param id: :class:`.Target` ID
        :returns: a :class:`.PoseSet` of the subset

        """
        assert isinstance(id, int)
        values = self.db.select_where(
            query="pose_id", table="pose", key="target", value=id, multiple=True
        )
        ids = [v for v, in values if v]

        target = self.db.get_target(id=id)

        pset = self[ids]
        pset._name = f'poses for "{target}"'
        return pset

    def get_by_smiles(self, smiles: str) -> "Pose | PoseSet | None":
        """Get a member pose by it's smiles"""

        from .tools import inchikey_from_smiles, sanitise_smiles, SanitisationError

        try:
            flat_smiles = sanitise_smiles(smiles, sanitisation_failed="error")
        except SanitisationError as e:
            mrich.error(f"Could not sanitise {smiles=}")
            mrich.error(str(e))
            return None
        except AssertionError:
            mrich.error(f"Could not sanitise {smiles=}")
            return None
            return c

        # get the compound

        flat_inchikey = inchikey_from_smiles(flat_smiles)

        comp_id = self.db.select_id_where(
            table="compound", key="inchikey", value=flat_inchikey
        )

        if not comp_id:
            return None

        (comp_id,) = comp_id

        # get the poses

        pose_ids = self.db.select_id_where(
            table="pose", key="compound", value=comp_id, multiple=True
        )

        if not pose_ids:
            return None

        pose_ids = [i for i, in pose_ids]
        pset = self[pose_ids]

        # identify the pose

        inchikey = inchikey_from_smiles(smiles)

        matches = set()
        for pose in pset:
            if pose.inchikey == inchikey:
                matches.add(pose.id)
        matches = list(matches)

        if not matches:
            mrich.error(f"Did not find pose matching stereochemistry (C{comp_id})")
            return None

        if len(matches) == 1:
            return self[matches[0]]

        return self[matches]

    def get_by_subsite(
        self,
        *,
        id: int,
    ) -> "PoseSet":
        """Get all child poses with a certain :class:`.Subsite` ID:

        :param id: :class:`.Subsite` ID
        :returns: a :class:`.PoseSet` of the subset

        """
        assert isinstance(id, int)
        values = self.db.select_where(
            query="subsite_tag_pose",
            table="subsite_tag",
            key="ref",
            value=id,
            multiple=True,
        )
        ids = [v for v, in values if v]

        subsite = self.db.get_subsite_name(id=id)

        pset = self[ids]
        pset._name = f'poses in "{subsite}"'
        return pset

    def get_by_metadata(
        self,
        key: str,
        value: str | None = None,
    ) -> "PoseSet":
        """Get all child poses by their metadata. If no value is passed, then simply containing the key in the metadata dictionary is sufficient

        :param key: metadata key to match
        :param value: metadata value to match, if ``None`` any pose with the key present will be returned (Default value = None)
        :returns: a :class:`.PoseSet` of the subset

        """
        results = self.db.select(
            query="pose_id, pose_metadata", table="pose", multiple=True
        )
        if value is None:
            ids = [i for i, d in results if d and f'"{key}":' in d]
            name = f"poses with {key} in metadata"
        else:
            if isinstance(value, str):
                value = f'"{value}"'
            ids = [i for i, d in results if d and f'"{key}": {value}' in d]
            name = f"poses with metadata[{key}] == {value}"

        pset = self[ids]
        pset._name = name
        return pset

    def get_by_metadata_substring_match(
        self,
        substring: str,
    ) -> "PoseSet":
        """Get :class:`.PoseSet` of poses with metadata JSON containing substring"""

        assert substring
        assert isinstance(substring, str)

        pose_ids = self.db.select_where(
            table="pose",
            query="pose_id",
            key=f"""pose_metadata LIKE '%{substring}%'""",
            multiple=True,
        )

        if not pose_ids:
            mrich.error("No poses with export ")
            return None

        pose_ids = [i for i, in pose_ids]

        name = f"poses with '{substring}' in metadata"

        pset = self[pose_ids]
        pset._name = name

        return pset

    def draw(
        self,
        max_draw: int = 100,
    ) -> None:
        """Render the poses

        :param max_draw: show a warning if trying to draw more than this number of poses (Default value = 100)

        """
        if len(self) <= max_draw:
            self[:].draw()
        else:
            mrich.warning(
                f"Too many poses: {len(self)} > {max_draw=}. Increase max_draw or use animal.poses[:].draw()"
            )

    def summary(self) -> None:
        """Print a summary of this pose set"""
        mrich.header("PoseTable()")
        mrich.var("#poses", len(self))
        mrich.var("tags", self.tags)

    def interactive(self) -> None:
        """Interactive widget to navigate poses in the table

        .. attention::

                This method instantiates a :class:`.PoseSet` containing all poses, it is recommended to instead select a subset for display. This method is only intended for use within a Jupyter Notebook.

        """

        self[self.ids].interactive()

    ### DUNDERS

    def __call__(
        self,
        *,
        tag: str | None = None,
        target: int | None = None,
        subsite: int | None = None,
        smiles: str | None = None,
    ) -> "PoseSet":
        """Filter poses by a given tag, subsite ID, or target ID. See :meth:`.PoseTable.get_by_tag`, :meth:`.PoseTable.get_by_target`, amd :meth:`.PoseTable.get_by_subsite`"""

        if tag:
            return self.get_by_tag(tag)
        elif target:
            return self.get_by_target(id=target)
        elif subsite:
            return self.get_by_subsite(id=subsite)
        elif smiles:
            return self.get_by_smiles(smiles=smiles)
        else:
            raise NotImplementedError

    def __getitem__(
        self,
        key: int | str | tuple | list | set | slice,
    ) -> Pose:
        """Get a member :class:`.Pose` object or subset :class:`.PoseSet` thereof.

        :param key: Can be an integer ID, negative integer index, alias or inchikey string, list/set/tuple of IDs, or slice of IDs

        """

        from pandas import Series

        match key:

            case int():
                if key == 0:
                    return self.__getitem__(key=1)

                if key < 0:
                    key = len(self) + 1 + key
                    return self.__getitem__(key=key)

                else:
                    return self.db.get_pose(id=key)

            case str():
                pose = self.db.get_pose(alias=key)
                if not pose:
                    pose = self.db.get_pose(inchikey=key)
                return pose

            case key if (
                isinstance(key, list)
                or isinstance(key, tuple)
                or isinstance(key, set)
                or isinstance(key, Series)
            ):

                indices = []
                for i in key:
                    if isinstance(i, int):
                        index = i
                    elif isinstance(i, str):
                        index = self.db.get_pose_id(alias=i)
                        if not index:
                            index = self.db.get_pose_id(inchikey=i)
                    else:
                        raise NotImplementedError

                    assert index
                    indices.append(index)

                return PoseSet(self.db, indices)

            case slice():
                ids, name = self.db.slice_ids(
                    table=self.table,
                    start=key.start,
                    stop=key.stop,
                    step=key.step,
                    name=True,
                )
                pset = self[ids]
                pset._name = name
                return pset

            case _:
                mrich.error(
                    f"Unsupported type for PoseTable.__getitem__(): {type(key)}"
                )

        return None

    def __str__(self):
        """Unformatted string representation"""
        if self.name:
            s = f"{self.name}: "
        else:
            s = ""

        s += "{" f"P × {len(self)}" "}"

        return s

    def __repr__(self) -> str:
        """ANSI Formatted string representation"""
        return f"{mcol.bold}{mcol.underline}{self}{mcol.unbold}{mcol.ununderline}"

    def __rich__(self) -> str:
        """Rich Formatted string representation"""
        return f"[bold underline]{self}"

    def __len__(self) -> int:
        """Total number of compounds"""
        return self.db.count(self.table)

    def __iter__(self):
        """Iterate through all compounds"""
        return iter(self[i + 1] for i in range(len(self)))


class PoseSet:
    """Object representing a subset of the 'pose' table in the :class:`.Database`.

    .. attention::

            :class:`.PoseSet` objects should not be created directly. Instead use the :meth:`.HIPPO.poses` property. See :doc:`getting_started` and :doc:`insert_elaborations`.

    Use as an iterable
    ==================

    Iterate through :class:`.Pose` objects in the set:

    ::

            pset = animal.poses[:100]

            for pose in pset:
                    ...

    Check membership
    ================

    To determine if a :class:`.Pose` is present in the set:

    ::

            is_member = pose in cset

    Selecting compounds in the set
    ==============================

    The :class:`.PoseSet` can be indexed like standard Python lists by their indices

    ::

            pset = animal.poses[1:100]

            # indexing individual compounds
            pose = pset[0]  # get the first pose
            pose = pset[1]  # get the second pose
            pose = pset[-1] # get the last pose

            # getting a subset of compounds using a slice
            pset2 = pset[13:18] # using a slice

    """

    _table = "pose"

    def __init__(
        self,
        db: Database,
        indices: list = None,
        *,
        sort: bool = True,
        name: str | None = None,
    ) -> None:

        self._db = db

        indices = indices or []

        if not isinstance(indices, list):
            indices = list(indices)

        assert all(isinstance(i, int) for i in indices)

        if sort:
            self._indices = sorted(list(set(indices)))
        else:

            # remove duplicates but keep order
            self._indices = dict()
            for i in indices:
                if i not in self._indices:
                    self._indices[i] = i
            self._indices = list(self._indices.keys())

        self._interactions = None

        self._name = name

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
    def indices(self) -> list[int]:
        """Returns the ids of poses in this set"""
        return self._indices

    @property
    def ids(self) -> list[int]:
        """Returns the ids of poses in this set"""
        return self._indices

    @property
    def name(self) -> str | None:
        """Returns the name of set"""
        return self._name

    @property
    def names(self) -> list[str]:
        """Returns the aliases of poses in this set"""
        return [p.name for p in self]

    @property
    def aliases(self) -> list[str]:
        """Returns the aliases of child poses"""
        return [
            self.db.select_where(
                table=self.table, query="pose_alias", key="id", value=i, multiple=False
            )[0]
            for i in self.indices
        ]

    @property
    def inchikeys(self) -> list[str]:
        """Returns the inchikeys of child poses"""
        return [
            self.db.select_where(
                table=self.table,
                query="pose_inchikey",
                key="id",
                value=i,
                multiple=False,
            )[0]
            for i in self.indices
        ]

    @property
    def id_name_dict(self) -> dict:
        """Return a dictionary mapping pose ID's to their name"""

        records = self.db.select_where(
            table=self.table,
            query="pose_id, pose_inchikey, pose_alias",
            key=f"pose_id IN {self.str_ids}",
            multiple=True,
        )

        lookup = {}
        for i, inchikey, alias in records:
            if alias:
                lookup[i] = alias
            else:
                lookup[i] = inchikey

        return lookup

    @property
    def smiles(self) -> list[str]:
        """Returns the smiles of poses in this set"""
        pairs = self.db.select_where(
            table=self.table,
            query="pose_id, pose_smiles",
            key=f"pose_id IN {self.str_ids}",
            multiple=True,
        )

        results = []
        for pose_id, smiles in pairs:
            if smiles is None:
                pose = self.db.get_pose(id=pose_id)
                smiles = pose.smiles

            results.append(smiles)

        return results

    @property
    def tags(self) -> set[str]:
        """Returns the set of unique tags present in this pose set"""
        values = self.db.select_where(
            table="tag",
            query="DISTINCT tag_name",
            key=f"tag_pose in {self.str_ids}",
            multiple=True,
        )
        return set(v for v, in values)

    @property
    def compounds(self) -> "CompoundSet":
        """Get the compounds associated to this set of poses"""
        from .cset import CompoundSet

        ids = self.db.select_where(
            table="pose",
            query="DISTINCT pose_compound",
            key=f"pose_id in {self.str_ids}",
            multiple=True,
        )
        ids = [v for v, in ids]
        return CompoundSet(self.db, ids)

    @property
    def mols(self) -> "list[rdkit.Chem.mol]":
        """Get the rdkit Molecules contained in this set"""
        return [p.mol for p in self]

    @property
    def num_compounds(self) -> int:
        """Count the compounds associated to this set of poses"""
        return len(self.compounds)

    @property
    def df(self) -> "pandas.DataFrame":
        """Get a DataFrame of the poses in this set"""
        return self.get_df(mol=True)

    @property
    def reference_ids(self) -> set[int]:
        """Return a set of :class:`.Pose` ID's of the all the distinct references in this :class:`.PoseSet`"""
        values = self.db.select_where(
            table="pose",
            query="DISTINCT pose_reference",
            key=f"pose_reference IS NOT NULL and pose_id in {self.str_ids}",
            value=None,
            multiple=True,
        )
        return set(v for v, in values)

    @property
    def inspiration_sets(self) -> list[set[int]]:
        """Return a list of unique sets of inspiration :class:`.Pose` IDs"""

        sql = f"""
        SELECT inspiration_derivative, inspiration_original FROM inspiration
        WHERE inspiration_derivative IN {self.str_ids}
        """

        pairs = self.db.execute(sql).fetchall()

        data = {}
        for derivative, original in pairs:
            if derivative not in data:
                data[derivative] = set()
            data[derivative].add(original)

        data = {k: tuple(sorted(list(v))) for k, v in data.items()}

        unique = set(data.values())

        return unique

    @property
    def num_inspiration_sets(self) -> int:
        """Return the number of unique sets of inspirations"""
        return len(self.inspiration_sets)

    @property
    def num_inspirations(self) -> int:
        """Return the number of unique inspirations for poses in this set"""
        (count,) = self.db.select_where(
            table="inspiration",
            query="COUNT(DISTINCT inspiration_original)",
            key=f"inspiration_derivative IN {self.str_ids}",
        )

        return count

    @property
    def inspirations(self) -> int:
        """Return the number of unique inspirations for poses in this set"""
        records = self.db.select_where(
            table="inspiration",
            query="DISTINCT inspiration_original",
            key=f"inspiration_derivative IN {self.str_ids}",
            multiple=True,
        )

        if not records:
            return None

        return PoseSet(self.db, [i for i, in records])

    @property
    def str_ids(self) -> str:
        """Return an SQL formatted tuple string of the :class:`.Pose` IDs"""
        return str(tuple(self.ids)).replace(",)", ")")

    @property
    def targets(self) -> "list[Target]":
        """Returns the :class:`.Target` objects of poses in this set"""
        return [self.db.get_target(id=q) for q in self.target_ids]

    @property
    def target_names(self) -> list[str]:
        """Returns the :class:`.Target` objects of poses in this set"""
        return [self.db.get_target_name(id=q) for q in self.target_ids]

    @property
    def target_ids(self) -> list[int]:
        """Returns the :class:`.Target` objects ID's of poses in this set"""
        result = self.db.select_where(
            table=self.table,
            query="DISTINCT pose_target",
            key=f"pose_id in {self.str_ids}",
            multiple=True,
        )
        return [q for q, in result]

    @property
    def best_placed_pose(self) -> Pose:
        """Returns the pose with the best distance_score in this subset"""
        return self.db.get_pose(id=self.best_placed_pose_id)

    @property
    def best_placed_pose_id(self) -> int:
        """Get the id of the pose with the best distance_score in this subset"""
        query = f"pose_id, MIN(pose_distance_score)"
        query = self.db.select_where(
            table="pose", query=query, key=f"pose_id in {self.str_ids}", multiple=False
        )
        return query[0]

    @property
    def interactions(self) -> "InteractionSet":
        """Get a :class:`.InteractionSet` for this :class:`.Pose`"""
        if not self._interactions:
            from .iset import InteractionSet

            self._interactions = InteractionSet.from_pose(self)
        return self._interactions

    @property
    def interaction_overlap_score(self) -> int:
        """Count the number of member pose pairs which share at least one but not all interactions"""

        sql = f"""
        SELECT DISTINCT interaction_pose, feature_id, interaction_type FROM interaction 
        INNER JOIN feature ON interaction_feature = feature_id
        WHERE interaction_pose IN {self.str_ids}
        """

        # mrich.print(sql)

        records = self.db.execute(sql).fetchall()

        ISETS = {}
        for pose_id, feature_id, interaction_type in records:
            values = ISETS.get(pose_id, set())
            values.add((interaction_type, feature_id))
            ISETS[pose_id] = values

        ids = [i for i in self.ids if i in ISETS]

        count = 0
        for pose_j in ids:
            iset_j = ISETS[pose_j]
            for pose_k in ids:
                iset_k = ISETS[pose_k]

                # try:
                # except KeyError:
                #     mrich.error("No interactions for pose with id", pose_k, "Has it been fingerprinted?")
                #     continue

                intersection = iset_j & iset_k
                diff1 = iset_j - iset_k
                diff2 = iset_k - iset_j

                if intersection and diff1 and diff2:
                    count += 1

        return count

    def get_interaction_clusters(self) -> "dict[int, PoseSet]":
        """Cluster poses based on shared interactions."""

        import networkx as nx
        import community as louvain
        from itertools import combinations

        # get interaction records

        sql = f"""
        SELECT DISTINCT interaction_pose, feature_residue_name, feature_residue_number, interaction_type FROM interaction 
        INNER JOIN feature ON interaction_feature = feature_id
        WHERE interaction_pose IN {self.str_ids}
        """

        records = self.db.execute(sql).fetchall()

        ISETS = {}
        for (
            pose_id,
            feature_residue_name,
            feature_residue_number,
            interaction_type,
        ) in records:
            values = ISETS.get(pose_id, set())
            values.add((interaction_type, feature_residue_name, feature_residue_number))
            ISETS[pose_id] = values

        pairs = combinations(ISETS.keys(), 2)

        # construct overlap dictionary

        OVERLAPS = {}
        for id1, id2 in pairs:
            iset1 = ISETS[id1]
            iset2 = ISETS[id2]
            OVERLAPS[(id1, id2)] = len(iset1 & iset2)

        # make the graph
        G = nx.Graph()

        for (id1, id2), count in OVERLAPS.items():
            G.add_edge(id1, id2, weight=count)

        # partition the graph

        partition = louvain.best_partition(G, weight="weight")

        # find the clusters

        clusters = {}
        for node, cluster_id in partition.items():
            clusters.setdefault(cluster_id, set()).add(node)

        # create the PoseSets

        psets = {
            i: PoseSet(self.db, ids, name=f"Cluster {i}")
            for i, ids in enumerate(clusters.values())
        }

        all_ids = set(sum((pset.ids for pset in psets.values()), []))

        # calculate modal interactions

        for i, cluster in psets.items():

            mrich.var(cluster.name, len(cluster), unit="poses")

            df = cluster.interactions.df

            unique_counts = df.groupby(["type", "residue_name", "residue_number"])[
                "pose_id"
            ].nunique()

            max_count = unique_counts.max()
            max_pairs = unique_counts[unique_counts == max_count]

            for (
                interaction_type,
                residue_name,
                residue_number,
            ) in max_pairs.index.values:
                mrich.print(interaction_type, "w/", residue_name, residue_number)

        # unclustered
        unclustered = set((i for i in self.ids if i not in all_ids))
        psets[None] = PoseSet(self.db, unclustered, name="Unclustered")

        return psets

    @property
    def num_fingerprinted(self) -> int:
        """Count the number of fingerprinted poses in this set"""
        return self.db.count_where(
            table="pose", key=f"pose_id IN {self.str_ids} AND pose_fingerprint = 1"
        )

    @property
    def fraction_fingerprinted(self) -> float:
        """Return the fraction of fingerprinted poses in this set"""
        return self.num_fingerprinted / len(self)

    @property
    def num_subsites(self) -> int:
        """Count the number of subsites that poses in this set come into contact with"""
        (count,) = self.db.select_where(
            query="COUNT(DISTINCT subsite_tag_ref)",
            table="subsite_tag",
            key=f"subsite_tag_pose IN {self.str_ids}",
            none="quiet",
        )
        if count is None:
            count = 0
        return count

    @property
    def subsite_balance(self) -> float:
        """Measure of how evenly subsite counts are distributed across poses in this set"""

        from numpy import std

        sql = f"""
        SELECT COUNT(DISTINCT subsite_tag_ref) FROM subsite_tag
        WHERE subsite_tag_pose IN {self.str_ids}
        GROUP BY subsite_tag_pose
        """

        counts = self.db.execute(sql).fetchall()

        counts = [c for c, in counts] + [0 for _ in range(len(self) - len(counts))]

        return -std(counts)

    @property
    def subsite_ids(self) -> set[int]:
        """Return a list of subsite id's of member poses"""

        sql = f"""
        SELECT DISTINCT subsite_tag_ref FROM subsite_tag
        WHERE subsite_tag_pose IN {self.str_ids}
        """

        subsite_ids = self.db.execute(sql).fetchall()

        if not subsite_ids:
            return set()

        subsite_ids = set([i for i, in subsite_ids])

        return subsite_ids

    @property
    def avg_energy_score(self) -> float:
        """Average energy score of poses in this set"""

        from numpy import mean

        sql = f"""
        SELECT pose_energy_score FROM pose
        WHERE pose_id IN {self.str_ids}
        """

        scores = self.db.execute(sql).fetchall()
        return mean([s for s, in scores if s is not None])

    @property
    def avg_distance_score(self) -> float:
        """Average distance score of poses in this set"""

        from numpy import mean

        sql = f"""
        SELECT pose_distance_score FROM pose
        WHERE pose_id IN {self.str_ids}
        """

        scores = self.db.execute(sql).fetchall()

        return mean([s for s, in scores if s is not None])

    @property
    def derivatives(self) -> "PoseSet":
        ids = self.db.select_where(
            table="inspiration",
            query="inspiration_derivative",
            key=f"inspiration_original IN {self.str_ids}",
            multiple=True,
            none="quiet",
        )
        if not ids:
            return None
        ids = [i for i, in ids]
        pset = PoseSet(self.db, ids, name=f"derivatives of {self}")
        return pset

    ### FILTERING

    def get_by_tag(
        self,
        tag: str,
        inverse: bool = False,
    ) -> "PoseSet":
        """Get all child poses with a certain tag

        :param tag: tag to filter by
        :param inverse: return all poses *not* tagged with ``tag`` (Default value = False)

        """
        values = self.db.select_where(
            query="tag_pose", table="tag", key="name", value=tag, multiple=True
        )
        if inverse:
            matches = [v for v, in values if v]
            ids = [i for i in self.ids if i not in matches]
        else:
            ids = [v for v, in values if v and v in self.ids]
        return PoseSet(self.db, ids)

    def get_by_metadata(
        self, key: str, value: str | None = None, debug: bool = False
    ) -> "PoseSet":
        """Get all child poses with by their metadata. If no value is passed, then simply containing the key in the metadata dictionary is sufficient

        :param key: metadata key to search for
        :param value: metadata value, if ``None`` return poses with the metadata key regardless of value (Default value = None)

        """
        results = self.db.select_where(
            query="pose_id, pose_metadata",
            key=f"pose_id IN {self.str_ids}",
            table="pose",
            multiple=True,
        )

        if value is None:
            ids = [i for i, d in results if d and f'"{key}":' in d]

        else:
            if isinstance(value, str):
                value = f'"{value}"'

            ids = []

            for i, d in results:
                if not d:
                    continue

                if debug:
                    mrich.print(i, d, f'"{key}": {value}' in d)

                if f'"{key}": {value}' in d:
                    ids.append(i)
                else:
                    continue

                if debug:
                    break

        return PoseSet(self.db, ids)

    def get_by_inspiration(self, inspiration: int | Pose, inverse: bool = False):
        """Get all child poses with with this inspiration.

        :param inspiration: inspiration :class:`.Pose` ID or object
        :param inverse: invert the selection (Default value = False)

        """

        ids = set()

        for pose in self:
            if not inverse:
                for pose_inspiration in pose.inspirations:
                    if pose_inspiration == inspiration:
                        ids.add(pose.id)
                        break

            elif inverse:
                for pose_inspiration in pose.inspirations:
                    if pose_inspiration == inspiration:
                        break
                else:
                    ids.add(pose.id)

        return PoseSet(self.db, ids)

    def get_df(
        self, skip_no_mol=True, reference: str = "name", mol: bool = False, **kwargs
    ) -> "pandas.DataFrame":
        """Get a DataFrame of the poses in this set. Keyword arguments passed to :meth:`.Pose.get_dict`.

        :param skip_no_mol: skip poses that have no mol (Default value = True)

        """

        from pandas import DataFrame

        data = []

        # if len(self) > 100:
        gen = mrich.track(enumerate(self), prefix="PoseSet --> DataFrame")
        track = True
        # else:
        #     gen = enumerate(self)
        #     track = False

        for i, pose in gen:

            d = pose.get_dict(reference=reference, mol=mol, **kwargs)

            mrich.set_progress_field("progress", f"{i+1}/{len(self)}")

            if skip_no_mol and not d["mol"]:
                mrich.warning(f'Skipping pose with no mol: {d["id"]} {d["name"]}')
                continue
            data.append(d)

        return DataFrame(data)

    def get_by_reference(
        self,
        ref_id: int,
    ) -> "PoseSet | None":
        """Get poses with a certain reference id

        :param ref_id: reference :class:`.Pose` ID

        """
        values = self.db.select_where(
            table="pose",
            query="pose_id",
            key=f"pose_reference={ref_id} AND pose_id in {self.str_ids}",
            multiple=True,
        )
        if not values:
            return None
        return PoseSet(self.db, [v for v, in values])

    def get_by_compound(
        self,
        *,
        compound: "int | Compound",
    ) -> "PoseSet | None":
        """Select a subset of this :class:`.PoseSet` by the associated :class:`.Compound`.

        :param compound: :class:`.Compound` object or ID
        :returns: a :class:`.PoseSet` of the selection

        """
        from .compound import Compound

        if isinstance(compound, Compound):
            compound = compound.id

        values = self.db.select_where(
            query="pose_id",
            table="pose",
            key=f"pose_compound={compound} AND pose_id in {self.str_ids}",
            multiple=True,
            none="quiet",
        )
        if not values:
            return None
        ids = [v for v, in values if v]
        return PoseSet(self.db, [v for v, in values])

    def get_by_target(
        self,
        *,
        id: int,
    ) -> "PoseSet | None":
        """Select a subset of this :class:`.PoseSet` by the associated :class:`.Target`.

        :param id: :class:`.Target` ID
        :returns: a :class:`.PoseSet` of the selection

        """
        assert isinstance(id, int)
        values = self.db.select_where(
            query="pose_id",
            table="pose",
            key=f"pose_target is {id} AND pose_id in {self.str_ids}",
            multiple=True,
            none="quiet",
        )
        ids = [v for v, in values if v]
        if not ids:
            return None
        return PoseSet(self.db, ids)

    def get_by_subsite(
        self,
        *,
        id: int,
    ) -> "PoseSet | None":
        """Select a subset of this :class:`.PoseSet` by the associated :class:`.Subsite`.

        :param id: :class:`.Subsite` ID
        :returns: a :class:`.PoseSet` of the selection

        """
        assert isinstance(id, int)
        values = self.db.select_where(
            query="subsite_tag_pose",
            table="subsite_tag",
            key=f"subsite_tag_ref is {id} AND subsite_tag_pose in {self.str_ids}",
            multiple=True,
            none="quiet",
        )
        ids = [v for v, in values if v]
        if not ids:
            return None

        if self.name:
            name = f"{self.name} & subsite={id}"
        else:
            name = None

        return PoseSet(self.db, ids, name=name)

    def filter(
        self,
        function,
        inverse: bool = False,
    ):
        """Filter this :class:`.PoseSet` by selecting members where ``function(pose)`` is truthy

        :param function: callable object
        :param inverse: invert the selection (Default value = False)

        """

        ids = set()
        for pose in self:
            value = function(pose)
            # mrich.debug(f'{pose=} {value=}')
            if value and not inverse:
                ids.add(pose.id)
            elif not value and inverse:
                ids.add(pose.id)

        return PoseSet(self.db, ids)

    ### BULK SETTING

    @property
    def reference(self):
        """Bulk set the references for poses in this set"""
        raise NotImplementedError(
            "This attribute only allows setting, ``PoseSet.reference = ...``"
        )

    @reference.setter
    def reference(self, r) -> None:
        """Bulk set the references for poses in this set"""
        if not isinstance(r, int):
            assert r._table == "pose"
            r = r.id

        for i in self.indices:
            self.db.update(
                table="pose", id=i, key="pose_reference", value=r, commit=False
            )

        self.db.commit()

    def add_tag(
        self,
        tag: str,
    ) -> None:
        """Add this tag to every member of the set"""

        assert isinstance(tag, str)

        for i in self.indices:
            self.db.insert_tag(name=tag, pose=i, commit=False)

        mrich.print(f'Tagged {self} w/ "{tag}"')

        self.db.commit()

    def append_to_metadata(
        self,
        key,
        value,
    ) -> None:
        """Append a specific item to list-like values associated with a given key for all member's metadata dictionaries

        :param key: the :class:`.Metadata` key to match
        :param value: the value to append to the list

        """
        for id in self.indices:
            metadata = self.db.get_metadata(table="pose", id=id)
            metadata.append(key, value)

    ### SPLITTING

    def split_by_reference(self) -> "dict[int,PoseSet]":
        """Split this :class:`.PoseSet` into subsets grouped by reference ID

        :returns: a dictionary with reference :class:`.Pose` IDs as keys and :class:`.PoseSet` subsets as values

        """
        sets = {}
        for ref_id in self.reference_ids:
            sets[ref_id] = self.get_by_reference(ref_id)
        return sets

    def split_by_inspirations(
        self,
        single_set: bool = False,
    ) -> "dict[int,PoseSet] | PoseSet":
        """Split this :class:`.PoseSet` into subsets grouped by inspirations

        :param single_set: Return a single :class:`.PoseSet` with members sorted by inspirations (Default value = False)
        :returns: a dictionary with tuples of inspiration :class:`.Pose` IDs as keys and :class:`.PoseSet` subsets as values

        """

        sets = {}

        for pose in self:

            insp_ids = tuple(pose.get_inspiration_ids())

            if insp_ids not in sets:
                sets[insp_ids] = PoseSet(self.db, [pose.id])
            else:
                sets[insp_ids]._indices.append(pose.id)

        if single_set:
            mrich.var("#unique inspiration combinations", len(sets))
            sets = PoseSet(self.db, sum([s.ids for s in sets.values()], []), sort=False)

        return sets

    ### EXPORTING

    def write_sdf(
        self,
        out_path: str,
        name_col: str = "name",
        inspirations: bool | str = "fragalysis",
        **kwargs,
    ) -> None:
        """Write an SDF

        :param out_path: filepath of the output
        :param name_col: pose property to use as the name column, can be ``["name", "alias", "inchikey", "id"]`` (Default value = 'name')
        :param inspirations: Include inspirations? ``[True, False, 'fragalysis']`` Specify ``fragalysis`` to format as a comma-separated string (Default value = "fragalysis")

        """

        from pathlib import Path
        import json

        df = self.get_df(mol=True, inspirations=inspirations, **kwargs)

        if name_col not in ["name", "alias", "inchikey", "id"]:
            # try getting name from metadata
            records = self.db.select_where(
                table="pose",
                query="pose_id, pose_metadata",
                key=f"pose_id IN {self.str_ids}",
                multiple=True,
            )

            longcode_lookup = {}
            for i, d in records:
                if d:
                    metadata = json.loads(d)
                else:
                    metadata = {}

                longcode_lookup[i] = metadata.get(name_col, None)

            values = []
            for i, row in df.iterrows():
                values.append(longcode_lookup[row["id"]])

            df[name_col] = values

        df.rename(inplace=True, columns={name_col: "_Name", "mol": "ROMol"})

        mrich.writing(out_path)

        from rdkit.Chem import PandasTools

        PandasTools.WriteSDF(df, out_path, "ROMol", "_Name", list(df.columns))

        # keep record of export
        value = str(Path(out_path).resolve())
        self.db.remove_metadata_list_item(table="pose", key="exports", value=value)
        self.append_to_metadata(key="exports", value=value)

    def to_fragalysis(
        self,
        out_path: str,
        *,
        method: str,
        ref_url: str = "https://hippo.winokan.com",
        submitter_name: str,
        submitter_email: str,
        submitter_institution: str,
        metadata: bool = True,
        sort_by: str | None = None,
        sort_reverse: bool = False,
        generate_pdbs: bool = False,
        ingredients: IngredientSet = None,
        skip_no_reference: bool = True,
        skip_no_inspirations: bool = True,
        skip_metadata: list[str] | None = None,
        tags: bool = True,
        extra_cols: dict[str, list] = None,
        name_col: str = "name",
        **kwargs,
    ):
        """Prepare an SDF for upload to the RHS of Fragalysis.

        :param out_path: the file path to write to
        :param method: method used to generate the compounds
        :param ref_url: reference URL for the method
        :param submitter_name: name of the person submitting the compounds
        :param submitter_email: email of the person submitting the compounds
        :param submitter_institution: institution name of the person submitting the compounds
        :param metadata: include metadata in the output? (Default value = True)
        :param skipmetadata: exclude metadata keys from output
        :param sort_by: if set will sort the SDF by this column/field (Default value = None)
        :param sort_reverse: reverse the sorting (Default value = False)
        :param generate_pdbs: generate accompanying protein-ligand complex PDBs (Default value = False)
        :param ingredients: get procurement and amount information from this :class:`.IngredientSet` (Default value = None)
        :param tags: include a column for tags in the output (Default value = True)
        :param extra_cols: extra_cols should be a dictionary with a key for each column name, and list values where the first element is the field description, and all subsequent elements are values for each pose.
        :param name: How to determine the molecule name, see :meth:`.PoseSet.get_df`

        """

        from .fragalysis import generate_header
        from pathlib import Path
        from rdkit.Chem import SDWriter, PandasTools

        assert out_path.endswith(".sdf")

        _name_col = "_Name"
        mol_col = "ROMol"

        # make sure references are defined:

        # values = self.db.select_where(
        #     table="pose",
        #     query="DISTINCT pose_id",
        #     key=f"pose_reference IS NULL and pose_id in {self.str_ids}",
        #     multiple=True,
        #     none="quiet",
        # )

        # if values:
        #     poses_missing_refs = set(v for v, in values)
        #     added_refs = PoseSet(self.db)

        #     for pose in PoseSet(self.db, poses_missing_refs):
        #         if "hits" in pose.tags:
        #             pose.reference = pose.id
        #             added_refs._indices.append(pose.id)

        #     poses_missing_refs -= set(added_refs.ids)
        # else:
        #     poses_missing_refs = None

        # if poses_missing_refs:
        #     mrich.warning(f"{len(poses_missing_refs)} Poses missing reference")
        #     mrich.var("poses w/o reference", poses_missing_refs)
        #     poses = PoseSet(self.db, set(self.ids) - poses_missing_refs)

        mrich.debug(len(self), "poses in set")
        poses = None

        if skip_no_reference:

            values = self.db.select_where(
                table="pose",
                query="DISTINCT pose_id",
                key=f"pose_reference IS NOT NULL and pose_id in {self.str_ids}",
                multiple=True,
                none="error",
            )

            poses = PoseSet(self.db, [i for i, in values])

            mrich.debug(len(poses), "remaining after skipping null reference")

        if skip_no_inspirations:

            if not poses:
                poses = self

            values = self.db.select_where(
                table="inspiration",
                query="DISTINCT inspiration_derivative",
                key=f"inspiration_derivative IN {poses.str_ids}",
                multiple=True,
                none="error",
            )

            poses = PoseSet(self.db, [i for i, in values])

            mrich.debug(len(poses), "remaining after skipping null inspirations")

        if not poses:
            poses = PoseSet(self.db, self.ids)

        mrich.var("#poses", len(poses))

        # get the dataframe of poses

        pose_df = poses.get_df(
            mol=True,
            inspirations="names",
            subsites="names",
            duplicate_name="original ID",
            reference="name",
            metadata=metadata,
            tags=tags,
            sanitise_null_metadata_values=True,
            sanitise_tag_list_separator=";",
            sanitise_metadata_list_separator=";",
            skip_metadata=skip_metadata,
            **kwargs,
        )

        drops = ["path", "compound", "target", "ref_pdb", "original SMILES"]

        if ingredients:
            drops.pop(drops.index("compound"))

        prev = len(pose_df)
        pose_df = pose_df[pose_df["reference"].notna()]
        if len(pose_df) < prev:
            mrich.warning(f"Skipping {prev - len(pose_df)} Poses with no reference")

        pose_df = pose_df.drop(columns=drops, errors="ignore")

        pose_df[_name_col] = pose_df[name_col]

        pose_df.rename(
            inplace=True,
            columns={
                "id": "HIPPO Pose ID",
                "compound_id": "HIPPO Compound ID",
                "mol": mol_col,
                "inspirations": "ref_mols",
                "reference": "ref_pdb",
                "smiles": "original SMILES",
                "compound": "compound inchikey",
            },
        )

        extras = {
            "HIPPO Pose ID": "HIPPO Pose ID",
            "HIPPO Compound ID": "HIPPO Compound ID",
            "smiles": "smiles",
            "ref_pdb": "protein reference",
            "ref_mols": "fragment inspirations",
            "original ID": "original ID",
            "compound inchikey": "compound inchikey",
            "distance_score": "distance_score",
            "energy_score": "energy_score",
            "subsites": "subsites",
        }

        if extra_cols:
            for key, value in extra_cols.items():
                extras[key] = value[0]

        if ingredients:

            q_entries = []
            q_prices = []
            q_lead_times = []
            q_amounts = []

            currency = None

            for i, row in pose_df.iterrows():

                compound_id = self.db.get_compound_id(inchikey=row["compound inchikey"])

                ingredient = ingredients(compound_id=compound_id)

                if isinstance(ingredient, IngredientSet):
                    ingredient = sorted(
                        [i for i in ingredient], key=lambda x: x.quote.price
                    )[0]

                quote = ingredient.quote
                if not currency:
                    currency = quote.currency
                else:
                    assert quote.currency == currency

                q_entries.append(quote.entry_str)
                q_prices.append(quote.price)
                q_lead_times.append(quote.lead_time)
                q_amounts.append(quote.amount)

            pose_df["Supplier Catalogue Entry"] = q_entries
            # pose_df['Supplier:Catalogue:Entry'] = q_entries
            pose_df[f"Price ({currency})"] = q_prices
            pose_df["Lead time (working days)"] = q_lead_times
            pose_df["Amount (mg)"] = q_amounts

            extras["Supplier Catalogue Entry"] = "Supplier Catalogue Entry string"
            extras[f"Price ({currency})"] = "Quoted price"
            extras["Lead time (working days)"] = "Quoted lead-time"
            extras["Amount (mg)"] = "Quoted amount"

        if generate_pdbs:

            from zipfile import ZipFile

            # output subdirectory
            out_key = Path(out_path).name.removesuffix(".sdf")
            pdb_dir = Path(out_path).parent / Path(out_key)
            pdb_dir.mkdir(exist_ok=True)
            zip_path = Path(out_path).parent / f"{out_key}_pdbs.zip"

            # create the zip archive
            with ZipFile(str(zip_path.resolve()), "w") as z:

                # loop over poses
                for (i, row), pose in zip(pose_df.iterrows(), poses):

                    # filenames
                    pdb_name = f"{out_key}_{row._Name}.pdb"
                    pdb_path = pdb_dir / pdb_name
                    pose_df.loc[i, "ref_pdb"] = pdb_name

                    # generate the PL-complex
                    sys = pose.complex_system

                    # write the PDB
                    mrich.writing(pdb_path)
                    sys.write(pdb_path, verbosity=0)
                    z.write(pdb_path)

            mrich.writing(f"{out_key}_pdbs.zip")

        # create the header molecule

        df_cols = set(pose_df.columns)

        header = generate_header(
            self[0],
            method=method,
            ref_url=ref_url,
            submitter_name=submitter_name,
            submitter_email=submitter_email,
            submitter_institution=submitter_institution,
            extras=extras,
            metadata=metadata,
        )

        header_cols = set(header.GetPropNames())

        # # empty properties
        # pose_df["generation_date"] = [None] * len(pose_df)
        # pose_df["submitter_name"] = [None] * len(pose_df)
        # pose_df["method"] = [None] * len(pose_df)
        # pose_df["submitter_email"] = [None] * len(pose_df)
        # pose_df["ref_url"] = [None] * len(pose_df)

        if extra_cols:
            for key, value in extra_cols.items():
                if len(value) != len(pose_df) + 1:
                    mrich.error(
                        f'extra_col "{key}" does not have the correct number of values'
                    )
                    raise ValueError(
                        f'extra_col "{key}" does not have the correct number of values'
                    )
                pose_df[key] = value[1:]

        if sort_by:
            pose_df = pose_df.sort_values(by=sort_by, ascending=not sort_reverse)

        fields = []

        mrich.writing(out_path)

        with open(out_path, "w") as sdfh:
            with SDWriter(sdfh) as w:
                w.write(header)
            PandasTools.WriteSDF(
                pose_df, sdfh, mol_col, _name_col, set(pose_df.columns)
            )

        # keep record of export
        value = str(Path(out_path).resolve())
        self.db.remove_metadata_list_item(table="pose", key="exports", value=value)
        self.append_to_metadata(key="exports", value=value)

        return pose_df

    def to_pymol(self, prefix: str | None = None) -> None:
        """Group the poses by reference protein and inspirations and output relevant PDBs and SDFs.

        :param prefix: prefix to give all output subdirectories (Default value = None)

        """

        commands = []

        prefix = prefix or ""
        if prefix:
            prefix = f"{prefix}_"

        from pathlib import Path

        for i, (ref_id, poses) in enumerate(self.split_by_reference().items()):

            ref_pose = self.db.get_pose(id=ref_id)
            ref_name = ref_pose.name or ref_id

            # create the subdirectory
            ref_dir = Path(f"{prefix}ref_{ref_name}")
            mrich.writing(ref_dir)
            ref_dir.mkdir(parents=True, exist_ok=True)

            # write the reference protein
            ref_pdb = ref_dir / f"ref_{ref_name}.pdb"
            ref_pose.protein_system.write(ref_pdb, verbosity=0)

            # color the reference:
            commands.append(f"load {ref_pdb.resolve()}")
            commands.append("hide")
            commands.append("show lines")
            commands.append("show surface")
            commands.append("util.cbaw")
            commands.append("set surface_color, white")
            commands.append("set transparency,  0.4")

            for j, (insp_ids, poses) in enumerate(
                poses.split_by_inspirations().items()
            ):

                inspirations = PoseSet(self.db, insp_ids)
                insp_names = "-".join(inspirations.names)

                # create the subdirectory
                insp_dir = ref_dir / insp_names
                insp_dir.mkdir(parents=True, exist_ok=True)

                # write the inspirations
                insp_sdf = insp_dir / f"{insp_names}_frags.sdf"
                inspirations.write_sdf(insp_sdf)

                commands.append(f"load {insp_sdf.resolve()}")
                commands.append(
                    f"set all_states, on, {insp_sdf.name.removesuffix('.sdf')}"
                )
                commands.append(
                    f"util.rainbow \"{insp_sdf.name.removesuffix('.sdf')}\""
                )

                # write the poses
                pose_sdf = insp_dir / f"{insp_names}_derivatives.sdf"
                poses.write_sdf(pose_sdf)

                commands.append(f"load {pose_sdf.resolve()}")
                commands.append(f'util.cbaw "{pose_sdf.name.removesuffix(".sdf")}"')

                if j > 0:
                    commands.append(f"disable \"{insp_sdf.name.removesuffix('.sdf')}\"")
                    commands.append(f'disable "{pose_sdf.name.removesuffix(".sdf")}"')

        return "; ".join(commands)

    def to_knitwork(
        self, out_path: str, path_root: str = ".", aligned_files_dir: str | None = None
    ) -> None:
        """Knitwork takes a CSV input with:

        - observation shortcode
        - smiles
        - path_to_ligand_mol
        - path_to_pdb

        :param out_path: path to output CSV
        :param path_root: paths in CSV will be relative to here

        """

        from os.path import relpath
        from pathlib import Path

        out_path = Path(out_path).resolve()
        path_root = Path(path_root).resolve()
        mrich.var("out_path", out_path)
        mrich.var("path_root", path_root)
        mrich.var("aligned_files_dir", aligned_files_dir)

        assert out_path.name.endswith(".csv")

        with open(out_path, "wt") as f:

            mrich.writing(out_path)

            for pose in self:

                assert pose.alias
                assert "hits" in pose.tags

                if aligned_files_dir:

                    mol = str(pose.mol_path)
                    pdb = str(pose.apo_path)

                    assert "aligned_files" in mol
                    assert "aligned_files" in pdb

                    mol = mol.split("aligned_files/")[-1]
                    pdb = pdb.split("aligned_files/")[-1]

                    aligned_files_dir = Path(aligned_files_dir)

                    mol = relpath(aligned_files_dir / mol, path_root)
                    pdb = relpath(aligned_files_dir / pdb, path_root)

                else:
                    mol = relpath(pose.mol_path, path_root)
                    pdb = relpath(pose.apo_path, path_root)

                data = [pose.alias, pose.compound.smiles, mol, pdb]

                f.write(",".join(data))
                f.write("\n")

    def to_syndirella(
        self, out_key: "str | Path", separate: bool = False
    ) -> "DataFrame":
        """Create syndirella inputs"""

        from pathlib import Path

        out_key = Path(out_key)

        if separate:
            dfs = []
            from pandas import concat

            for i, pose in enumerate(self):
                mrich.h3(f"{i}/{len(self)}: {pose}")
                this_out_key = out_key / f"P{pose.id}"
                df = pose.to_syndirella(out_key=this_out_key)
                dfs.append(df)
            return concat(dfs)

        import shutil
        from pandas import DataFrame

        out_dir = out_key.parent
        out_key = out_key.name

        out_dir.mkdir(parents=True, exist_ok=True)

        template_dir = out_dir / "templates"
        mrich.writing(template_dir)
        template_dir.mkdir(parents=True, exist_ok=True)

        # mrich.var("out_dir", out_dir)
        mrich.var("out_key", out_key)
        mrich.var("#poses", len(self))

        pset_inspirations = set()

        data = []
        for pose in mrich.track(self, prefix="preparing inputs..."):

            comp = pose.compound

            ref = pose.reference
            inspirations = pose.inspirations

            if not ref:
                mrich.warning(pose, "has no reference, using self as template")
                ref = pose
                assert ref.apo_path, f"Reference {ref} has no apo_path"

            if not inspirations:
                mrich.warning(pose, "has no inspirations, using self")
                inspirations = PoseSet(self.db, [pose.id])

            for i in inspirations.ids:
                pset_inspirations.add(i)

            d = dict(
                smiles=comp.smiles,
                # template=ref.path,
                template=ref.name,
                compound_set=out_key,
            )

            mrich.writing(template_dir / ref.apo_path.name)
            shutil.copy(ref.apo_path, template_dir / ref.apo_path.name)

            for i, p in enumerate(inspirations):
                d[f"hit{i+1}"] = p.name

            data.append(d)

        df = DataFrame(data)

        csv_name = out_dir / f"{out_key}_syndirella_input.csv"
        mrich.writing(csv_name)
        df.to_csv(csv_name, index=False)

        inspirations = PoseSet(self.db, pset_inspirations)

        sdf_name = out_dir / f"{out_key}_syndirella_inspiration_hits.sdf"
        inspirations.write_sdf(
            sdf_name,
            inspirations=False,
            tags=False,
            metadata=False,
            reference=False,
            name_col="name",
        )

        return df

    ### OUTPUT

    def interactive(
        self,
        print_name: str = True,
        method: str | None = None,
        function: Callable | None = None,
        **kwargs,
    ):
        """Interactive widget to navigate compounds in the table

        :param print_name: print the :class:`.Pose` name  (Default value = True)
        :param method: pass the name of a :class:`.Pose` method to interactively display. Keyword arguments to interactive() will be passed through (Default value = None)
        :param function: pass a callable which will be called as `function(pose)`

        """

        from ipywidgets import (
            interactive,
            BoundedIntText,
            Checkbox,
            interactive_output,
            HBox,
            GridBox,
            Layout,
            VBox,
        )
        from IPython.display import display
        from pprint import pprint

        if method:

            def widget(i):
                pose = self[i]
                if print_name:
                    print(repr(pose))
                value = getattr(pose, method)(**kwargs)
                if value:
                    display(value)

            return interactive(
                widget,
                i=BoundedIntText(
                    value=0,
                    min=0,
                    max=len(self) - 1,
                    step=1,
                    description="Pose:",
                    disabled=False,
                ),
            )

        elif function:

            def widget(i):
                pose = self[i]
                if print_name:
                    display(pose)
                function(pose)

            return interactive(
                widget,
                i=BoundedIntText(
                    value=0,
                    min=0,
                    max=len(self) - 1,
                    step=1,
                    description="Pose:",
                    disabled=False,
                ),
            )

        else:

            a = BoundedIntText(
                value=0,
                min=0,
                max=len(self) - 1,
                step=1,
                description=f"Pose (/{len(self)}):",
                disabled=False,
            )

            b = Checkbox(description="Name", value=True)
            c = Checkbox(description="Summary", value=False)
            h = Checkbox(description="Tags", value=False)
            i = Checkbox(description="Subsites", value=False)
            d = Checkbox(description="2D (Comp.)", value=False)
            e = Checkbox(description="2D (Pose)", value=False)
            f = Checkbox(description="3D", value=True)
            g = Checkbox(description="Metadata", value=False)

            ui1 = GridBox(
                [b, c, d, h],
                layout=Layout(grid_template_columns="repeat(4, 100px)"),
            )
            ui2 = GridBox(
                [e, f, g, i],
                layout=Layout(grid_template_columns="repeat(4, 100px)"),
            )
            ui = VBox([a, ui1, ui2])

            def widget(
                i,
                name=True,
                summary=True,
                grid=True,
                draw2d=True,
                draw=True,
                tags=True,
                subsites=True,
                metadata=True,
            ):
                pose = self[i]
                if name:
                    print(repr(pose))

                if summary:
                    pose.summary(metadata=False, tags=False, subsites=False)
                if tags:
                    print(pose.tags)
                if subsites:
                    print(pose.subsites)
                if grid:
                    pose.grid()
                if draw2d:
                    pose.draw2d()
                if draw:
                    pose.draw()
                if metadata:
                    mrich.title("Metadata:")
                    pprint(pose.metadata)

            out = interactive_output(
                widget,
                {
                    "i": a,
                    "name": b,
                    "summary": c,
                    "grid": d,
                    "draw2d": e,
                    "draw": f,
                    "metadata": g,
                    "tags": h,
                    "subsites": i,
                },
            )

            display(ui, out)

    def summary(self) -> None:
        """Print a summary of this pose set"""
        mrich.header("PoseSet()")
        mrich.var("#poses", len(self))
        mrich.var("#compounds", self.num_compounds)
        mrich.var("tags", self.tags)

    def draw(self) -> None:
        """Render this pose set with Py3Dmol"""

        from molparse.rdkit import draw_mols

        mols = [p.mol for p in self]

        drawing = draw_mols(mols)
        # display(drawing)

    def grid(self) -> None:
        """Draw a grid of all contained molecules"""
        from molparse.rdkit import draw_grid

        data = [(p.name, p.compound.mol) for p in self]

        mols = [d[1] for d in data]
        labels = [d[0] for d in data]

        drawing = draw_grid(mols, labels=labels)
        display(drawing)

    ### PRIVATE

    def _delete(self, *, force: bool = False) -> None:
        """Delete poses in this set"""

        if not force:
            mrich.warning("Deleting Poses is risky! Set force=True to continue")
            return

        # delete the poses in this set
        self.db.delete_where(table=self.table, key=f"pose_id IN {self.str_ids}")

        # check for other references to this pose
        self.db.delete_where(table="tag", key=f"tag_pose IN {self.str_ids}")
        self.db.delete_where(
            table="inspiration", key=f"inspiration_original IN {self.str_ids}"
        )
        self.db.delete_where(
            table="inspiration", key=f"inspiration_derivative IN {self.str_ids}"
        )
        self.db.delete_where(
            table="subsite_tag", key=f"subsite_tag_pose IN {self.str_ids}"
        )
        self.db.delete_where(
            table="interaction", key=f"interaction_pose IN {self.str_ids}"
        )

        self.db.execute(
            f"""
            UPDATE pose
            SET pose_reference = NULL
            WHERE pose_id IN {self.str_ids}
        """
        )

    ### DUNDERS

    def __str__(self):
        """Unformatted string representation"""
        if self.name:
            s = f"{self.name}: "
        else:
            s = ""

        s += "{" f"P × {len(self)}" "}"

        return s

    def __repr__(self) -> str:
        """ANSI Formatted string representation"""
        return f"{mcol.bold}{mcol.underline}{self}{mcol.unbold}{mcol.ununderline}"

    def __rich__(self) -> str:
        """Rich Formatted string representation"""
        return f"[bold underline]{self}"

    def __len__(self) -> int:
        """The number of poses in this set"""
        return len(self.indices)

    def __iter__(self):
        """Iterate through poses in this set"""
        return iter(self.db.get_pose(id=i) for i in self.indices)

    def __getitem__(
        self,
        key: int | slice,
    ) -> "Pose | PoseSet":
        """Get poses or subsets thereof from this set

        :param key: integer index or slice of indices

        """
        match key:

            case int():
                try:
                    index = self.indices[key]
                except IndexError:
                    mrich.error(f"list index out of range: {key=} for {self}")
                    raise
                return self.db.get_pose(id=index)

            case slice():
                ids = self.indices[key]
                return PoseSet(self.db, ids)

            case _:
                raise NotImplementedError

    def __add__(
        self,
        other: "PoseSet",
    ) -> "PoseSet":
        """Add a :class:`.PoseSet` to this set"""
        if isinstance(other, PoseSet):
            return PoseSet(self.db, self.ids + other.ids, sort=False)
        elif isinstance(other, Pose):
            return PoseSet(self.db, self.ids + [other.id], sort=False)
        else:
            raise NotImplementedError

    def __sub__(
        self,
        other: "PoseSet",
    ) -> "PoseSet":
        """Substract a :class:`.PoseSet` from this set"""
        match other:
            case PoseSet():
                ids = set(self.ids) - set(other.ids)
                return PoseSet(self.db, ids, sort=False)
            case int():
                # assert other in set(self.ids)
                return PoseSet(self.db, [i for i in self.ids if i != other], sort=False)

    def __call__(
        self,
        *,
        tag: str = None,
        target: int = None,
        subsite: int = None,
    ) -> "PoseSet":
        """Filter poses by a given tag, Subsite ID, or target ID. See :meth:`.PoseSet.get_by_tag`, :meth:`.PoseSet.get_by_target`, amd :meth:`.PoseSet.get_by_subsite`"""

        if tag:
            return self.get_by_tag(tag)
        elif target:
            return self.get_by_target(id=target)
        elif subsite:
            return self.get_by_subsite(id=subsite)
        else:
            raise NotImplementedError
