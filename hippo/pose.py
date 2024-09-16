import mcol
import molparse as mp

from rdkit import Chem

import numpy as np

from .tags import TagSet

import pickle

from pathlib import Path

import logging

logger = logging.getLogger("HIPPO")

from molparse.rdkit.features import (
    FEATURE_FAMILIES,
    COMPLEMENTARY_FEATURES,
    INTERACTION_TYPES,
)

INTERACTION_CUTOFF = {
    "Hydrophobic": 4.5,
    "Hydrogen Bond": 3.5,
    "Electrostatic": 4.5,
    "π-stacking": 6.0,
    "π-cation": 4.5,
}

PI_STACK_MIN_CUTOFF = 3.8
PI_STACK_F2F_CUTOFF = 4.5
PI_STACK_E2F_CUTOFF = 6.0


class Pose:
    """A :class:`.Pose` is a particular conformer of a :class:`.Compound` within a protein environment. A pose will have its own (stereochemical) smiles string, and must have a path to a coordinate file. Poses can have *inspirations* that can be used to trace fragment-derived scaffolds in merges and expansions.

    .. attention::

            :class:`.Pose` objects should not be created directly. Instead use :meth:`.HIPPO.register_pose` or :meth:`.HIPPO.poses`

    """

    _table = "pose"

    def __init__(
        self,
        db: "Database",
        id: int,
        inchikey: str | None,
        alias: str | None,
        smiles: str,
        reference: int,  # another pose
        path: str,
        compound: int,
        target: int,
        mol: Chem.Mol | bytes | None,
        fingerprint: int,
        energy_score: float | None = None,
        distance_score: float | None = None,
        metadata: dict | None = None,
    ):

        self._db = db
        self._id = id
        self._inchikey = inchikey
        self._alias = alias
        self._smiles = smiles
        self._compound_id = compound
        self._target = target
        self._path = path
        self._protein_system = None
        self._energy_score = energy_score
        self._distance_score = distance_score

        self._base_id = None
        self._num_heavy_atoms = None

        self._has_fingerprint = False

        if fingerprint is None:
            self._has_fingerprint = False
        elif not isinstance(fingerprint, int):
            logger.warning("Legacy fingerprint data format")
            self.has_fingerprint = False
        else:
            self.has_fingerprint = bool(fingerprint)

        # print(f'{self}{metadata=}')
        self._metadata = metadata
        self._tags = None
        self._reference = reference
        self._reference_id = reference
        self._interactions = None

        if isinstance(mol, bytes):
            self._mol = Chem.Mol(mol)
        else:
            self._mol = mol

    ### FACTORIES

    ### PROPERTIES

    @property
    def db(self) -> "Database":
        """Returns a pointer to the parent database"""
        return self._db

    @property
    def id(self) -> int:
        """Returns the pose's database ID"""
        return self._id

    @property
    def inchikey(self) -> str:
        """Returns the pose's inchikey"""
        if not self._inchikey:
            self.smiles
        return self._inchikey

    @property
    def alias(self) -> str:
        """Returns the pose's alias"""
        return self._alias

    @property
    def name(self) -> str:
        """Returns the pose's name"""
        if n := self.alias:
            return n
        else:
            return self.inchikey

    @alias.setter
    def alias(self, n) -> None:
        """Set the pose's alias"""
        assert isinstance(n, str)
        self._alias = n
        self.db.update(table="pose", id=self.id, key="pose_alias", value=n)

    @inchikey.setter
    def inchikey(self, n) -> None:
        """Set the pose's inchikey"""
        assert isinstance(n, str)
        self._inchikey = n
        self.db.update(table="pose", id=self.id, key="pose_inchikey", value=n)

    @property
    def smiles(self) -> str:
        """Returns the pose's smiles"""
        if not self._smiles:
            from molparse.rdkit import mol_to_smiles
            from rdkit.Chem.inchi import MolToInchiKey

            try:
                mol = self.mol
                self._smiles = mol_to_smiles(mol)
                self.inchikey = MolToInchiKey(mol)
                self.db.update(
                    table="pose", id=self.id, key="pose_smiles", value=self._smiles
                )
            except InvalidMolError:
                logger.warning(f"Taking smiles from {self.compound}")
                self._smiles = self.compound.smiles
        return self._smiles

    @property
    def target(self) -> "Target":
        """Returns the pose's associated target"""
        if isinstance(self._target, int):
            self._target = self.db.get_target(id=self._target)
        return self._target

    @property
    def compound_id(self) -> int:
        """Returns the pose's associated compound ID"""
        return self._compound_id

    @property
    def compound(self) -> "Compound":
        """Returns the pose's associated compound"""
        return self.get_compound()

    @property
    def path(self) -> str:
        """Returns the pose's path"""
        return self._path

    @property
    def reference(self) -> "Pose":
        """Returns the pose's protein reference (another pose)"""
        if isinstance(self._reference, int):
            self._reference = self.db.get_pose(id=self._reference)
        return self._reference

    @property
    def reference_id(self) -> int:
        """Returns the pose's protein reference ID"""
        return self._reference_id

    @reference.setter
    def reference(self, p):
        """Set the pose's reference"""
        if not isinstance(p, int):
            assert p._table == "pose"
            p = p.id
        self._reference = p
        self._reference_id = p
        self.db.update(table="pose", id=self.id, key="pose_reference", value=p)

    @property
    def mol(self) -> "rdkit.Chem.Mol":
        """Returns a pose's rdkit.Chem.Mol"""
        if not self._mol and self.path:

            if self.path.endswith(".pdb"):

                logger.reading(self.path)

                # logger.reading(self.path)
                sys = mp.parse(self.path, verbosity=False)

                self.protein_system = sys.protein_system

                # look for ligand mol from Fragalysis
                mol_path = list(
                    Path(self.path).parent.glob("*_ligand.mol")
                )  # str(Path(self.path).name).replace('.pdb','_ligand.mol')

                if len(mol_path) == 1:
                    mol_path = mol_path[0]
                    from rdkit.Chem import MolFromMolFile

                    mol = MolFromMolFile(str(mol_path.resolve()))

                elif len(mol_path) == 0:

                    lig_residues = sys["rLIG"]

                    if not lig_residues:
                        lig_residues = [r for r in sys.residues if r.type == "LIG"]

                    if len(lig_residues) > 1:
                        logger.warning(f"Multiple ligands in PDB {self}")

                    lig_res = lig_residues[0]

                    if not (mol := lig_res.rdkit_mol):
                        logger.error(
                            f"[{self}] Error computing RDKit Mol from PDB={self.path}"
                        )

                        print(lig_res.pdb_block)

                        lig_res.plot3d()

                        raise InvalidMolError

                    # clean up bond orders
                    from rdkit.Chem.AllChem import (
                        MolFromSmiles,
                        AssignBondOrdersFromTemplate,
                    )

                    template = MolFromSmiles(self.compound.smiles)
                    try:
                        mol = AssignBondOrdersFromTemplate(template, mol)
                    except Exception as e:
                        logger.error(
                            f"Exception occured during AssignBondOrdersFromTemplate for {self}.mol"
                        )
                        print(f"template_smiles={self.compound.smiles}")
                        print(f"pdbblock={print(lig_res.pdb_block)}")
                        logger.error(e)
                        mol = lig_res.rdkit_mol

                else:

                    logger.warning(
                        f"There are multiple *_ligand.mol files in {Path(self.path).parent}"
                    )

                self.mol = mol

            elif self.path.endswith(".mol"):

                logger.reading(self.path)

                # logger.reading(self.path)
                mol = mp.parse(self.path, verbosity=False)

                if not mol:
                    logger.error(
                        f"[{self}] Error computing RDKit Mol from .mol={self.path}"
                    )

                    raise InvalidMolError

                self.mol = mol

            else:

                raise NotImplementedError

            if not mol:
                logger.error(f"Could not parse {self}.path={self.path}")

        return self._mol

    @mol.setter
    def mol(self, m):
        """Set the pose's rdkit.Chem.Mol"""
        assert m
        from .tools import sanitise_mol

        self._mol = sanitise_mol(m)
        self.db.update(table="pose", id=self.id, key="pose_mol", value=m.ToBinary())

    @property
    def protein_system(self) -> "molparse.System":
        """Returns the pose's protein molparse.System"""
        if self._protein_system is None and self.path.endswith(".pdb"):
            # logger.debug(f'getting pose protein system {self}')
            self.protein_system = mp.parse(self.path, verbosity=False).protein_system
        return self._protein_system

    @protein_system.setter
    def protein_system(self, a):
        """Sets the pose's protein molparse.System"""
        self._protein_system = a

    @property
    def complex_system(self) -> "molparse.System":

        if self.has_complex_pdb_path:
            return mp.parse(self.path, verbosity=False)

        elif self.reference:

            # construct from .mol and reference

            system = self.reference.protein_system.copy()

            system.name = (
                f"{self.target.name}_{self.reference.name}_{self.compound.name}"
            )

            from molparse.rdkit import mol_to_AtomGroup

            ligand = mol_to_AtomGroup(self.mol)

            for atom in ligand.atoms:
                system.add_atom(atom)

            return system

        else:

            raise NotImplementedError

    @property
    def has_complex_pdb_path(self) -> bool:
        """Does this pose have a PDB file?"""
        return self.path.endswith(".pdb")

    @property
    def metadata(self) -> "MetaData":
        """Returns the pose's metadata"""
        if self._metadata is None:
            self._metadata = self.db.get_metadata(table="pose", id=self.id)
        return self._metadata

    @property
    def has_fingerprint(self) -> bool:
        """Does the pose have a fingerprint?"""
        return self._has_fingerprint

    @has_fingerprint.setter
    def has_fingerprint(self, fp):

        assert isinstance(fp, bool)

        self._has_fingerprint = fp

        # store in the database
        self.db.update(table="pose", id=self.id, key=f"pose_fingerprint", value=int(fp))

    @property
    def tags(self) -> "TagSet":
        """Returns the pose's tags"""
        if not self._tags:
            self._tags = self.get_tags()
        return self._tags

    @property
    def inspirations(self) -> "PoseSet":
        """Returns the pose's inspirations"""
        return self.get_inspirations()

    @property
    def features(self) -> "list[molparse.rdkit.Feature]":
        """Returns the pose's features"""
        return mp.rdkit.features_from_mol(self.mol)

    @property
    def dict(self) -> dict:
        """Serialised dictionary representing the pose"""
        return self.get_dict()

    @property
    def table(self) -> str:
        """Get the name of the database table"""
        return self._table

    @property
    def num_heavy_atoms(self) -> int:
        """Number of heavy atoms"""
        if not self._num_heavy_atoms:
            self._num_heavy_atoms = self.db.get_compound_computed_property(
                "num_heavy_atoms", self.compound_id
            )
        return self._num_heavy_atoms

    @property
    def num_atoms_added(self) -> int:
        """Calculate the number of atoms added relative to the base or inspirations"""
        if b_id := self.base_id:
            return self.num_atoms_added_wrt_base
        else:
            return self.num_atoms_added_wrt_inspirations

    @property
    def num_atoms_added_wrt_base(self) -> int | None:
        """Calculate the number of atoms added relative to the base"""
        if not (b_id := self.base_id):
            logger.error(f"{self} has no base")
            return None

        n_e = self.num_heavy_atoms
        n_b = self.db.get_compound_computed_property("num_heavy_atoms", b_id)
        return n_e - n_b

    @property
    def num_atoms_added_wrt_inspirations(self) -> int | None:
        """Calculate the number of atoms added relative to its inspirations"""

        inspirations = self.inspirations

        if not inspirations:
            logger.error(f"{self} has no inspirations")
            return None

        count = 0
        for i in inspirations:
            count += i.num_heavy_atoms

        if (self_count := self.num_heavy_atoms) is None:
            return None

        return self.num_heavy_atoms - count

    @property
    def base_id(self) -> int:
        """Get the base compound ID"""
        if not self._base_id:
            (val,) = self.db.select_where(
                table="compound",
                query="compound_base",
                key="id",
                value=self.compound_id,
            )
            self._base_id = val
        return self._base_id

    # @property
    # def fields(self):
    #   """ """
    #   return [p for p in dir(self) if not p.startswith('_')]

    @property
    def energy_score(self) -> float | None:
        """Energy score of the Pose (kcal/mol)"""
        return self._energy_score

    @property
    def distance_score(self) -> float | None:
        """Distance score of the Pose (w.r.t. its inspirations), in Angstroms"""
        return self._distance_score

    @property
    def interactions(self) -> "InteractionSet":
        """Get a :class:`.InteractionSet` for this :class:`.Pose`"""
        if not self._interactions:
            from .iset import InteractionSet

            self._interactions = InteractionSet.from_pose(self)
        return self._interactions

    @property
    def classic_fingerprint(self) -> dict:
        """Classic HIPPO fingerprint dictionary, mapping protein :class:`.Feature` ID's to the number of corresponding ligand features (from any :class:`.Pose`)"""
        return self.interactions.classic_fingerprint

    @property
    def subsites(self):

        from .subsite import SubsiteTag

        records = self.db.select_where(
            table="subsite_tag",
            key="pose",
            value=self.id,
            multiple=True,
            query="subsite_tag_id, subsite_tag_ref",
        )

        subsite_tags = []
        for record in records:
            id, ref = record
            subsite_tag = SubsiteTag(db=self.db, id=id, subsite_id=ref, pose_id=self.id)
            subsite_tags.append(subsite_tag)

        return subsite_tags

    ### METHODS

    def score_inspiration(
        self,
        debug: bool = False,
        draw: bool = False,
    ) -> float:
        """Score how well this Pose recapitulates the pharmacophoric features of its inspirations.

        :param debug: Increased verbosity for debugging (Default value = False)
        :param draw: Render each inspiration pose with it's features, the derivative with the combined features of the inspirations, and the derivative with it's features. (Default value = False)

        """

        from molparse.rdkit import SuCOS_score

        multi_sucos = SuCOS_score(
            self.inspirations.mols, self.mol, print_scores=debug, draw=draw
        )

        if debug:
            logger.var("energy_score", self.energy_score)
            logger.var("distance_score", self.distance_score)

            for inspiration in self.inspirations:
                logger.var(
                    f"{inspiration} SuCOS",
                    SuCOS_score(inspiration.mol, self.mol, print_scores=debug),
                )

            logger.var(f"multi SuCOS", multi_sucos)

        return multi_sucos

    def get_compound(self) -> "Compound":
        """Get the :class:`.Compound` that this pose is a conformer of"""
        return self.db.get_compound(id=self._compound_id)

    def get_tags(self) -> "TagSet":
        """Get this Pose's tags"""
        tags = self.db.select_where(
            query="tag_name",
            table="tag",
            key="pose",
            value=self.id,
            multiple=True,
            none="quiet",
        )
        return TagSet(self, {t[0] for t in tags})

    def get_inspiration_ids(self) -> list[int]:
        """Get the :class:`.Pose` IDs of this pose's inspirations"""
        inspirations = self.db.select_where(
            query="inspiration_original",
            table="inspiration",
            key="derivative",
            value=self.id,
            multiple=True,
            none="quiet",
        )
        if not inspirations:
            return None
        return set([v for v, in inspirations])

    def get_inspirations(self) -> "PoseSet":
        """Get a :class:`.PoseSet` of this pose's inspirations"""
        if not (inspirations := self.get_inspiration_ids()):
            return None

        from .pset import PoseSet

        return PoseSet(self.db, indices=inspirations)

    def get_dict(
        self,
        mol: bool = False,
        inspirations: bool | str = True,
        reference: bool | str = True,
        metadata: bool = True,
        duplicate_name: str | bool = False,
        tags: bool = True,
    ) -> dict:
        """Returns a dictionary representing this Pose. Arguments:

        :param mol: Include a ``rdkit.Chem.Mol`` in the output?  (Default value = False)
        :param inspirations: Include inspirations? ``[True, False, 'fragalysis']`` Specify ``fragalysis`` to format as a comma-separated string (Default value = True)
        :param reference: Include reference? ``[True, False, 'name']`` Specify ``name`` to include the :class:`.Pose` name rather than it's ID (Default value = True)
        :param metadata: Include metadata? (Default value = True)
        :param duplicate_name: Specify the name of a new column duplicating the pose name column  (Default value = False)
        :param tags: bool: Include tags? (Default value = True)

        """

        serialisable_fields = ["id", "inchikey", "alias", "name", "smiles", "path"]

        data = {}
        for key in serialisable_fields:
            data[key] = getattr(self, key)

        if duplicate_name:
            assert isinstance(duplicate_name, str)
            data[duplicate_name] = data["name"]

        if mol:
            try:
                data["mol"] = self.mol
            except InvalidMolError:
                data["mol"] = None

        data["compound"] = self.compound.name
        data["target"] = self.target.name

        if tags:
            data["tags"] = self.tags

        if inspirations == "fragalysis":
            if not self.inspirations:
                data["inspirations"] = ""
            else:
                data["inspirations"] = ",".join([p.name for p in self.inspirations])
        elif inspirations:
            data["inspirations"] = self.inspirations

        if reference == "name":
            if not self.reference:
                data["reference"] = ""
            else:
                data["reference"] = self.reference.name
        elif reference:
            data["reference"] = self.reference

        if metadata and (metadict := self.metadata):
            for key in metadict:
                data[key] = metadict[key]

        return data

    def add_subsite(self, name: str, commit: bool = True) -> "SubsiteTag":
        """Tag this pose with a protein subsite

        :param name: the name of the subsite
        :param commit: commit the insertion to the database
        :returns: :class:`.SubsiteTag`

        """

        id = self.db.insert_subsite_tag(pose_id=self.id, name=name, commit=commit)

        if not id:
            return None

        return self.db.get_subsite_tag(id=id)

    def calculate_interactions(
        self,
        resolve: bool = True,
        distance_padding: float = 0.0,
        angle_padding: float = 0.0,
        force: bool = False,
        debug: bool = False,
        commit: bool = True,
    ) -> None:
        """Enumerate all valid interactions between this ligand and the protein

        :param resolve: Cull duplicate / less-significant interactions
        :param distance_padding: Apply a padding in Angstrom to all distance cutoffs
        :param angle_padding: Apply a padding in degrees to all angle cutoffs
        :param force: Force a recalculation even if the pose has already been fingerprinted
        :param debug: Increase verbosity for debugging
        :param commit: commit the changes to the database (Default value = True)

        """

        if not self.has_fingerprint or force:

            def norm(coords):
                import numpy as np

                coords = np.array(coords).T
                cov = np.cov(coords)
                eig = np.linalg.eig(cov)
                vec = eig[1][:, 0]
                return vec

            def unit_vector(vector):
                """Returns the unit vector of the vector."""
                return vector / np.linalg.norm(vector)

            def angle_between(v1, v2):
                """Returns the angle in radians between vectors 'v1' and 'v2'::

                >>> angle_between((1, 0, 0), (0, 1, 0))
                1.5707963267948966
                >>> angle_between((1, 0, 0), (1, 0, 0))
                0.0
                >>> angle_between((1, 0, 0), (-1, 0, 0))
                3.141592653589793
                """
                v1_u = unit_vector(v1)
                v2_u = unit_vector(v2)
                a = 180 * np.arccos(np.clip(np.dot(v1_u, v2_u), -1.0, 1.0)) / np.pi

                if a > 90:
                    a = 180 - a
                return a

            # calculate interactions...

            ### clear old interactions

            self.db.delete_where(
                table="interaction", key="pose", value=self.id, commit=commit
            )
            self.has_fingerprint = False

            ### create temporary table

            self.db.create_table_interaction(table="temp_interaction", debug=False)

            ### load the ligand structure

            if self.path.endswith(".pdb"):
                from molparse import parse

                protein_system = self.protein_system
                if not self.protein_system:
                    protein_system = parse(self.path, verbosity=False).protein_system

            elif self.path.endswith(".mol") and self.reference:
                protein_system = self.reference.protein_system

            else:
                logger.debug("Unsupported: Pose.calculate_interactions()")
                raise NotImplementedError(f"{self.reference=}, {self.path=}")

            assert protein_system

            if not self.mol:
                logger.error(f"Could not read molecule for {self}")
                return

            ### get features

            comp_features = self.features
            protein_features = self.target.calculate_features(
                protein_system, reference_id=self.reference_id
            )

            ### organise ligand features by family
            comp_features_by_family = {}
            for family in FEATURE_FAMILIES:
                comp_features_by_family[family] = [
                    f for f in comp_features if f.family == family
                ]

            ### protein chain names
            chains = protein_system.chain_names

            # loop over protein features
            for prot_feature in protein_features:

                # skip chains that aren't present
                if prot_feature.chain_name not in chains:
                    continue

                prot_family = prot_feature.family

                prot_residue = protein_system.get_chain(
                    prot_feature.chain_name
                ).residues[f"n{prot_feature.residue_number}"]

                if not prot_residue:
                    continue

                if prot_residue.name != prot_feature.residue_name:
                    logger.warning(f"Feature {repr(prot_feature)}")
                    continue

                ### calculate protein coordinate
                prot_atoms = [
                    prot_residue.get_atom(a) for a in prot_feature.atom_names.split(" ")
                ]

                prot_coords = [a.np_pos for a in prot_atoms if a is not None]

                prot_coord = np.array(np.sum(prot_coords, axis=0) / len(prot_atoms))

                complementary_families = COMPLEMENTARY_FEATURES[prot_family]

                # print(prot_family, complementary_families)

                for complementary_family in complementary_families:

                    interaction_type = INTERACTION_TYPES[
                        (prot_family, complementary_family)
                    ]

                    complementary_comp_features = comp_features_by_family[
                        complementary_family
                    ]

                    for lig_feature in complementary_comp_features:

                        distance = np.linalg.norm(lig_feature - prot_coord)
                        angle = None

                        # check distance cutoff
                        if (
                            distance
                            > INTERACTION_CUTOFF[interaction_type] + distance_padding
                        ):
                            continue

                        # special rules for aromatics
                        if interaction_type.startswith("π"):
                            lig_coords = [
                                self.mol.GetConformer().GetAtomPosition(i - 1)
                                for i in lig_feature.atom_numbers
                            ]

                        # special rules for pi-stacking
                        if interaction_type == "π-stacking":

                            # calculate minimum distance
                            min_distance = None
                            for lig_coord in lig_coords:
                                for p_coord in prot_coords:
                                    d = np.linalg.norm(lig_coord - p_coord)
                                    if not min_distance or d < min_distance:
                                        min_distance = d

                            # skip interaction if no atom is within PI_STACK_MIN_CUTOFF
                            if min_distance > PI_STACK_MIN_CUTOFF + distance_padding:
                                if debug:
                                    print(
                                        f"skipping {prot_feature} due to pi-stack min_distance"
                                    )
                                    print(
                                        prot_feature.residue_name,
                                        prot_feature.residue_number,
                                        prot_feature.chain_name,
                                        prot_feature.family,
                                        lig_feature,
                                        distance,
                                    )
                                continue

                            ### angles
                            lig_norm = norm([list(p) for p in lig_coords])
                            prot_norm = norm(prot_coords)
                            angle = angle_between(lig_norm, prot_norm)

                            # Face to Face has more stringent restraints
                            if (
                                angle < 40 - angle_padding
                                and distance > PI_STACK_F2F_CUTOFF + distance_padding
                            ):
                                if debug:
                                    print(prot_feature.res_name_number_family_str)
                                    print(angle, distance)
                                continue

                        # special rules for pi-cation
                        elif interaction_type == "π-cation":

                            # construct vectors
                            if prot_family == "Aromatic":
                                aromatic_norm = norm(prot_coords)
                                cation_vec = lig_feature.position - prot_coord
                            else:
                                aromatic_norm = norm([list(p) for p in lig_coords])
                                cation_vec = prot_coord - lig_feature.position

                            # calculate angle
                            angle = angle_between(aromatic_norm, cation_vec)

                            # skip if angle too large
                            if angle > 30 + angle_padding:
                                if debug:
                                    print(prot_feature.res_name_number_family_str)
                                    print(angle, distance)
                                continue

                        if debug:
                            print(prot_feature, lig_feature)

                        # insert into the Database
                        self.db.insert_interaction(
                            feature=prot_feature.id,
                            pose=self.id,
                            type=interaction_type,
                            family=lig_feature.family,
                            atom_ids=lig_feature.atom_numbers,
                            prot_coord=prot_coord,
                            lig_coord=lig_feature.position,
                            distance=distance,
                            angle=angle,
                            energy=None,
                            commit=commit,
                            table="temp_interaction",
                        )

            self.has_fingerprint = True

            if resolve:
                from .iset import InteractionSet

                interactions = InteractionSet.from_pose(self, table="temp_interaction")
                interactions.resolve(debug=debug)
                # self.interactions.resolve(debug=debug, table='temp_interaction')

            ### transfer interactions from temporary table
            self.db.copy_temp_interactions()

            ### delete temporary table

            self.db.execute("DROP TABLE temp_interaction")

        elif debug:
            logger.warning(f"{self} is already fingerprinted, no new calculation")

    def calculate_classic_fingerprint(
        self,
        debug: bool = False,
    ) -> dict:
        """Calculate the pose's interaction fingerprint"""

        if self.path.endswith(".pdb"):

            import molparse as mp

            protein_system = self.protein_system
            if not self.protein_system:
                # logger.reading(self.path)
                protein_system = mp.parse(self.path, verbosity=False).protein_system

        elif self.path.endswith(".mol") and self.reference:

            # logger.debug('fingerprint from .mol and reference pose')
            protein_system = self.reference.protein_system

        else:

            logger.debug("Unsupported: Pose.calculate_fingerprint()")
            raise NotImplementedError(f"{self.reference=}, {self.path=}")

        assert protein_system

        if not self.mol:
            return

        comp_features = self.features

        comp_features_by_family = {}
        for family in FEATURE_FAMILIES:
            comp_features_by_family[family] = [
                f for f in comp_features if f.family == family
            ]

        # protein_features = self.target.features
        # if not protein_features:
        protein_features = self.target.calculate_features(protein_system)

        fingerprint = {}

        chains = protein_system.chain_names

        for prot_feature in protein_features:

            if prot_feature.chain_name not in chains:
                continue

            prot_family = prot_feature.family

            prot_residue = protein_system.get_chain(prot_feature.chain_name).residues[
                f"n{prot_feature.residue_number}"
            ]

            if not prot_residue:
                continue

            # if prot_feature.residue_number == 77:
            #   logger.debug(repr(prot_feature))

            if prot_residue.name != prot_feature.residue_name:
                logger.warning(f"Feature {repr(prot_feature)}")
                continue

            prot_atoms = [
                prot_residue.get_atom(a) for a in prot_feature.atom_names.split(" ")
            ]

            prot_coords = [a.np_pos for a in prot_atoms if a is not None]

            prot_coord = np.array(np.sum(prot_coords, axis=0) / len(prot_atoms))

            complementary_family = COMPLEMENTARY_FEATURES[prot_family]

            complementary_comp_features = comp_features_by_family[complementary_family]

            cutoff = FEATURE_PAIR_CUTOFFS[f"{prot_family} {complementary_family}"]

            valid_features = [
                f
                for f in complementary_comp_features
                if np.linalg.norm(f - prot_coord) <= cutoff
            ]

            if valid_features:
                if debug:
                    logger.debug(
                        f"PROT: {prot_feature.residue_name} {prot_feature.residue_number} {prot_feature.atom_names}, LIG: #{len(valid_features)} {[f for f in valid_features]}"
                    )
                fingerprint[prot_feature.id] = len(valid_features)

        return fingerprint

        # self.fingerprint = fingerprint

    def draw(
        self,
        inspirations: bool = True,
        protein: bool = False,
        **kwargs,
    ) -> None:
        """Render this pose (and its inspirations)

        :param inspirations: Render the inspirations? (Default value = True)
        :param protein: Render the protein? This wraps :meth:`.Pose.render` (Default value = False)

        """

        if protein:
            self.render(**kwargs)

        from molparse.rdkit import draw_mols

        mols = [self.mol]
        if inspirations and self.inspirations:
            mols += [i.mol for i in self.inspirations]

        draw_mols(mols)

    def render(
        self,
        protein="cartoon",
        ligand="stick",
        protein_color="spectrum",
    ) -> None:
        """Render this pose with the protein using py3Dmol

        :param protein: protein representation, default = 'cartoon'
        :param ligand: ligand representation, default = 'stick'
        :param protein_color: color of protein representation, default = 'spectrum'

        """

        from molparse.py3d import render

        display(
            render(
                self.complex_system,
                protein=protein,
                ligand=ligand,
                protein_color=protein_color,
            )
        )

    def grid(self) -> None:
        """Draw a grid of this pose with its inspirations"""
        from molparse.rdkit import draw_grid
        from IPython.display import display

        mols = [self.compound.mol]
        labels = [self.plain_repr()]
        if self.inspirations:
            mols += [i.compound.mol for i in self.inspirations]
            labels += [i.plain_repr() for i in self.inspirations]

        display(draw_grid(mols, labels=labels))

    def summary(self, metadata: bool = True) -> None:
        """Print a summary of this pose

        :param metadata: include metadata (Default value = True)

        """
        if self.alias:
            logger.header(f"{str(self)}: {self.alias}")
        else:
            logger.header(f"{str(self)}: {self.inchikey}")
        logger.var("inchikey", self.inchikey)
        logger.var("alias", self.alias)
        logger.var("smiles", self.smiles)
        logger.var("compound", repr(self.compound))
        logger.var("path", self.path)
        logger.var("target", repr(self.target))
        logger.var("reference", self.reference)
        logger.var("tags", self.tags)
        logger.var("num_heavy_atoms", self.num_heavy_atoms)
        if inspirations := self.inspirations:
            logger.var("inspirations", self.inspirations.names)
            logger.var("num_atoms_added", self.num_atoms_added)
        if metadata:
            logger.var("metadata", str(self.metadata))

    def showcase(self) -> None:
        """Print and render this pose as if you were using :meth:`.PoseSet.interactive`"""

        self.summary(metadata=False)
        self.grid()
        self.draw()
        from pprint import pprint

        logger.title("Metadata:")
        pprint(self.metadata)

    def plain_repr(self) -> str:
        """Unformatted detailed string representation"""
        if self.name:
            return f"{self.compound}->{self}: {self.name}"
        else:
            return f"{self.compound}->{self}"

    def plot3d(
        self,
        features: bool = False,
        **kwargs,
    ) -> "plotly.graph_objects.Figure":
        """Use Molparse/Plotly to create a 3d figure of this pose

        :param features: include the features in the figure
        :returns: a plotly Figure object

        """

        mol = self.mol

        import molparse as mp

        group = mp.rdkit.mol_to_AtomGroup(mol)

        if features:
            features = self.features

        return mp.go.plot3d(atoms=group.atoms, features=features, **kwargs)

    ### DUNDERS

    def __str__(self) -> str:
        """Unformatted string representation"""
        return f"P{self.id}"

    def __repr__(self) -> str:
        """Formatted string representation"""
        return f"{mcol.bold}{mcol.underline}{self.plain_repr()}{mcol.unbold}{mcol.ununderline}"

    def __eq__(self, other: "Pose") -> bool:
        """Compare this pose with another instance"""

        if isinstance(other, int):
            return self.id == other

        return self.id == other.id


class InvalidMolError(Exception):
    """Exception to be thrown when the molecule could not be parsed"""

    ...
