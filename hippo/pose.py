import mcol
import molparse as mp

from rdkit import Chem

import numpy as np

from .tags import TagSet

import pickle

from pathlib import Path

import mrich
from mrich import print

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
    "Sulfur-Sulfur": 4.0,  # https://pubs.acs.org/doi/full/10.1021/acs.cgd.5b01058
}

PI_STACK_MIN_CUTOFF = 3.8
PI_STACK_F2F_CUTOFF = 4.5
PI_STACK_E2F_CUTOFF = 6.0
MUTATION_WARNING_DIST = 15


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
        inspiration_score: float | None = None,
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
        self._inspiration_score = inspiration_score

        self._scaffold_ids = None
        self._num_heavy_atoms = None

        self._has_fingerprint = False

        if fingerprint is None:
            self._has_fingerprint = False
        elif not isinstance(fingerprint, int):
            mrich.warning("Legacy fingerprint data format")
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

        self._total_changes = db.total_changes
        self._num_atoms_added_wrt_inspirations = None

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
                mrich.warning(f"Taking smiles from {self.compound}")
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

                mrich.reading(self.path)

                # mrich.reading(self.path)
                sys = mp.parse(self.path, verbosity=False)

                self.protein_system = sys.protein_system

                sdf_path = list(Path(self.path).parent.glob("*_ligand.sdf"))

                if len(sdf_path) == 1:

                    supplier = Chem.SDMolSupplier(sdf_path[0])
                    mols = [mol for mol in supplier if mol is not None]

                    if len(mols) > 1:
                        mrich.warning(f"Multiple molecules in SDF {self}")

                    self.mol = mols[0]
                    return self._mol

                # look for ligand mol from Fragalysis
                mol_path = list(
                    Path(self.path).parent.glob("*_ligand.mol")
                )  # str(Path(self.path).name).replace('.pdb','_ligand.mol')

                if len(mol_path) == 1:
                    mol_path = mol_path[0]
                    from rdkit.Chem import MolFromMolFile

                    self._mol_path = mol_path.resolve()

                    mol = MolFromMolFile(str(self._mol_path))

                elif len(mol_path) == 0:

                    lig_residues = sys["rLIG"]

                    if not lig_residues:
                        lig_residues = [r for r in sys.residues if r.type == "LIG"]

                    if len(lig_residues) > 1:
                        mrich.warning(f"Multiple ligands in PDB {self}")

                    lig_res = lig_residues[0]

                    if not (mol := lig_res.rdkit_mol):
                        mrich.error(
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
                        mrich.error(
                            f"Exception occured during AssignBondOrdersFromTemplate for {self}.mol"
                        )
                        print(f"template_smiles={self.compound.smiles}")
                        print(f"pdbblock={print(lig_res.pdb_block)}")
                        mrich.error(e)
                        mol = lig_res.rdkit_mol

                else:

                    mrich.warning(
                        f"There are multiple *_ligand.mol files in {Path(self.path).parent}"
                    )

                self.mol = mol

            elif self.path.endswith(".mol"):

                mrich.reading(self.path)

                # mrich.reading(self.path)
                mol = mp.parse(self.path, verbosity=False)

                if not mol:
                    mrich.error(
                        f"[{self}] Error computing RDKit Mol from .mol={self.path}"
                    )

                    raise InvalidMolError

                self.mol = mol

            else:

                raise NotImplementedError

            if not mol:
                mrich.error(f"Could not parse {self}.path={self.path}")

        return self._mol

    @mol.setter
    def mol(self, m):
        """Set the pose's rdkit.Chem.Mol"""
        assert m
        from .tools import sanitise_mol

        self._mol = sanitise_mol(m)
        self.db.update(table="pose", id=self.id, key="pose_mol", value=m.ToBinary())

    @property
    def protonated_mol(self) -> "rdkit.Chem.Mol":
        """Guess hydrogen positions"""
        from rdkit.Chem import AllChem

        mol = self.mol
        protonated_mol = Chem.AddHs(mol)
        try:
            protonated_mol = AllChem.ConstrainedEmbed(protonated_mol, mol)
        except Exception as e:
            mrich.error("Error while embedding protonated molecule")
            mrich.error(e)
            return mol
        return protonated_mol

    @property
    def protein_system(self) -> "molparse.System":
        """Returns the pose's protein molparse.System"""
        if self._protein_system is None and self.path.endswith(".pdb"):
            # mrich.debug(f'getting pose protein system {self}')
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
        self.set_has_fingerprint(fp)

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
    def derivatives(self) -> "PoseSet":
        """Returns the pose's derivatives"""
        return self.get_derivatives()

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
        """Calculate the number of atoms added relative to the scaffold or inspirations"""
        if self.num_scaffolds == 1:
            return self.num_atoms_added_wrt_scaffolds
        else:
            return self.num_atoms_added_wrt_inspirations

    @property
    def num_atoms_added_wrt_scaffolds(self) -> int | list[int] | None:
        """Calculate the number of atoms added relative to the scaffold"""
        return self.compound.num_atoms_added

    @property
    def num_atoms_added_wrt_inspirations(self) -> int | None:
        """Calculate the number of atoms added relative to its inspirations"""

        if self._num_atoms_added_wrt_inspirations is None or self._db_changed:

            sql = f"""
            WITH inspirations AS (
                SELECT SUM(mol_num_hvyatms(compound_mol)) AS sum, inspiration_derivative FROM inspiration
                INNER JOIN pose ON inspiration_original = pose_id
                INNER JOIN compound ON pose_compound = compound_id
                WHERE inspiration_derivative = {self.id}
            )
            SELECT mol_num_hvyatms(compound_mol) - sum FROM inspirations
            INNER JOIN pose ON inspiration_derivative = pose_id
            INNER JOIN compound ON compound_id = pose_compound
            """

            (diff,) = self.db.execute(sql).fetchone()

            self._num_atoms_added_wrt_inspirations = diff

        return self._num_atoms_added_wrt_inspirations

    @property
    def num_scaffolds(self) -> int:
        """Get the number of scaffold scaffolds"""
        return len(self.scaffold_ids)

    @property
    def scaffold_ids(self) -> list[int] | None:
        """Get the scaffold :class:`.Compound` IDs"""
        if self._scaffold_ids is None:
            records = self.db.select_where(
                table="scaffold",
                query="scaffold_base",
                key="superstructure",
                value=self.compound_id,
                multiple=True,
                none="quiet",
            )
            records = [i for i, in records]
            self._scaffold_ids = records
        return self._scaffold_ids

    @property
    def energy_score(self) -> float | None:
        """Energy score of the Pose (kcal/mol)"""
        return self._energy_score

    @property
    def distance_score(self) -> float | None:
        """Distance score of the Pose (w.r.t. its inspirations), in Angstroms"""
        return self._distance_score

    @property
    def inspiration_score(self) -> float | None:
        """inspiration score of the Pose in range 0.00-1.00"""
        return self._inspiration_score

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
            none="quiet",
        )

        if not records:
            return None

        subsite_tags = []
        for record in records:
            id, ref = record
            subsite_tag = SubsiteTag(db=self.db, id=id, subsite_id=ref, pose_id=self.id)
            subsite_tags.append(subsite_tag)

        return subsite_tags

    @property
    def _db_changed(self) -> bool:
        """Has the database changed?"""
        if self._total_changes != self.db.total_changes:
            self._total_changes = self.db.total_changes
            return True
        return False

    @property
    def mol_path(self):
        path = Path(self.path)
        if path.name.endswith(".pdb"):
            mol_path = path.parent / path.name.replace("_hippo.pdb", ".pdb").replace(
                ".pdb", "_ligand.mol"
            )
            if not mol_path.exists():
                mol_path = path.parent / path.name.replace(
                    "_hippo.pdb", ".pdb"
                ).replace(".pdb", "_ligand.sdf")
                if not mol_path.exists():
                    mrich.error("Could not find ligand mol/sdf:", mol_path)
                    return None
            return mol_path
        elif path.name.endswith(".mol"):
            return path
        else:
            raise NotImplementedError

    @property
    def apo_path(self):
        path = Path(self.path)
        if path.name.endswith(".pdb"):
            apo_path = path.parent / path.name.replace("_hippo.pdb", ".pdb").replace(
                ".pdb", "_apo-desolv.pdb"
            )
            if not apo_path.exists():
                return None
            return apo_path
        else:
            raise NotImplementedError

    ### METHODS

    def score_inspiration(
        self,
        debug: bool = False,
        draw: bool = False,
        return_all: bool = False,
    ) -> float:
        """Score how well this Pose recapitulates the pharmacophoric features of its inspirations.

        :param debug: Increased verbosity for debugging (Default value = False)
        :param draw: Render each inspiration pose with it's features, the derivative with the combined features of the inspirations, and the derivative with it's features. (Default value = False)

        """

        # from molparse.rdkit import SuCOS_score
        from mucos import MuCOS_score

        multi_sucos = MuCOS_score(
            self.inspirations.mols,
            self.mol,
            print_scores=debug,
            draw=draw,
            return_all=return_all,
        )

        if debug:
            mrich.var("energy_score", self.energy_score)
            mrich.var("distance_score", self.distance_score)

            for inspiration in self.inspirations:
                mrich.var(
                    f"{inspiration} SuCOS",
                    MuCOS_score(inspiration.mol, self.mol, print_scores=debug),
                )

            mrich.var(f"multi SuCOS", multi_sucos)

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

    def get_derivative_ids(self) -> list[int]:
        """Get the :class:`.Pose` IDs of this pose's derivatives"""
        derivatives = self.db.select_where(
            query="inspiration_derivative",
            table="inspiration",
            key="original",
            value=self.id,
            multiple=True,
            none="quiet",
        )
        if not derivatives:
            return None
        return set([v for v, in derivatives])

    def get_inspirations(self) -> "PoseSet":
        """Get a :class:`.PoseSet` of this pose's inspirations"""
        if not (inspirations := self.get_inspiration_ids()):
            return None

        from .pset import PoseSet

        return PoseSet(self.db, indices=inspirations)

    def get_derivatives(self) -> "PoseSet":
        """Get a :class:`.PoseSet` of this pose's derivatives"""
        if not (derivatives := self.get_derivative_ids()):
            return None

        from .pset import PoseSet

        return PoseSet(self.db, indices=derivatives)

    def get_dict(
        self,
        mol: bool = False,
        inspirations: bool | str = True,
        subsites: bool | str = True,
        reference: bool | str = True,
        metadata: bool = True,
        duplicate_name: str | bool = False,
        sanitise_null_metadata_values: bool = False,
        skip_metadata: list[str] | None = None,
        sanitise_tag_list_separator: str | None = None,
        sanitise_metadata_list_separator: str | None = ";",
        tags: bool = True,
    ) -> dict:
        """Returns a dictionary representing this Pose. Arguments:

        :param mol: Include a ``rdkit.Chem.Mol`` in the output?  (Default value = False)
        :param inspirations: Include inspirations? ``[True, False, 'names']`` Specify ``names`` to format as a comma-separated string (Default value = True)
        :param subsites: Include subsites? ``[True, False, 'names']`` Specify ``names`` to format as a comma-separated string (Default value = True)
        :param reference: Include reference? ``[True, False, 'name']`` Specify ``name`` to include the :class:`.Pose` name rather than it's ID (Default value = True)
        :param metadata: Include metadata? (Default value = True)
        :param duplicate_name: Specify the name of a new column duplicating the pose name column  (Default value = False)
        :param tags: bool: Include tags? (Default value = True)

        """

        skip_metadata = skip_metadata or []

        serialisable_fields = [
            "id",
            "inchikey",
            "alias",
            "name",
            "smiles",
            "path",
            "distance_score",
            "energy_score",
            "inspiration_score",
        ]

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
        data["compound_id"] = self.compound.id
        data["target"] = self.target.name

        if tags:
            data["tags"] = self.tags
            if sanitise_tag_list_separator:
                data["tags"] = sanitise_tag_list_separator.join(data["tags"])

        if inspirations == "names":
            if not self.inspirations:
                data["inspirations"] = None
            else:
                data["inspirations"] = ",".join([p.name for p in self.inspirations])
        elif inspirations:
            data["inspirations"] = self.inspirations

        if subsites == "names":
            if not (sites := self.subsites):
                data["subsites"] = None
            else:
                data["subsites"] = ",".join([p.name for p in sites])
        elif subsites:
            data["subsites"] = self.subsites

        if reference == "name":
            if not self.reference:
                data["reference"] = ""
            else:
                data["reference"] = self.reference.name
        elif reference:
            data["reference"] = self.reference

        if metadata and (metadict := self.metadata):
            for key in metadict:

                value = metadict[key]

                if key in skip_metadata:
                    continue

                if (
                    sanitise_null_metadata_values
                    and isinstance(value, str)
                    and not value
                ):
                    value = None

                elif sanitise_metadata_list_separator and isinstance(value, list):

                    new_values = []

                    for v in value:
                        if (
                            sanitise_null_metadata_values
                            and isinstance(v, str)
                            and not v
                        ):
                            v = None

                        else:
                            v = str(v)

                        new_values.append(v)

                    value = sanitise_metadata_list_separator.join(new_values)

                data[key] = value

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
        mutation_warnings: bool = True,
        delete_temp_table: bool = True,
    ) -> None:
        """Enumerate all valid interactions between this ligand and the protein

        :param resolve: Cull duplicate / less-significant interactions
        :param distance_padding: Apply a padding in Angstrom to all distance cutoffs
        :param angle_padding: Apply a padding in degrees to all angle cutoffs
        :param force: Force a recalculation even if the pose has already been fingerprinted
        :param debug: Increase verbosity for debugging
        :param commit: commit the changes to the database (Default value = True)
        :param mutation_warnings: warn when there has been a mutation in the protein (Default value = True)
        :param delete_temp_table: delete the temporary interaction table created during interaction resolution (Default value = True)

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
            self.set_has_fingerprint(False, commit=commit)
            self._interactions = None

            ### create temporary table

            if "temp_interaction" in self.db.table_names:
                self.db.execute("DROP TABLE temp_interaction")

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
                mrich.debug("Unsupported: Pose.calculate_interactions()")
                raise NotImplementedError(f"{self}, {self.reference=}, {self.path=}")

            assert protein_system

            if not self.mol:
                mrich.error(f"Could not read molecule for {self}")
                return

            ### get features

            comp_features = self.features
            protein_features = self.target.calculate_features(
                protein_system, reference_id=self.reference_id
            )

            if debug:
                print("ligand features", comp_features)

            ### organise ligand features by family
            comp_features_by_family = {}
            for family in FEATURE_FAMILIES:
                comp_features_by_family[family] = [
                    f for f in comp_features if f.family == family
                ]

            ### protein chain names
            chains = protein_system.chain_names

            mutation_warnings = set()
            mutation_count = 0

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
                    com = prot_residue.centre_of_mass()
                    if any(
                        np.linalg.norm(com - cf.position) < MUTATION_WARNING_DIST
                        for cf in comp_features
                    ):
                        mutation_warnings.add(
                            f"{prot_residue.name} {prot_residue.number} -> {prot_feature.residue_name} {prot_feature.residue_number}"
                        )
                        mutation_count += 1
                    continue

                ### calculate protein coordinate
                prot_atoms = []
                for atom_name in prot_feature.atom_names.split(" "):
                    atom = prot_residue.get_atom(atom_name, verbosity=0)
                    if atom:
                        prot_atoms.append(atom)

                prot_coords = [a.np_pos for a in prot_atoms if a is not None]

                if not prot_coords:
                    # mrich.warning("Skipping feature with no atoms")
                    continue

                prot_coord = np.array(np.sum(prot_coords, axis=0) / len(prot_atoms))

                if prot_family not in COMPLEMENTARY_FEATURES:
                    continue

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
                            print("Prot:", prot_feature, "Lig:", lig_feature)

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

            if mutation_warnings:
                mrich.warning(
                    f"Skipped {mutation_count} protein features because the residue was mutated:"
                )
                for mutation in mutation_warnings:
                    mrich.warning(mutation)

            if resolve:
                from .iset import InteractionSet

                interactions = InteractionSet.from_pose(self, table="temp_interaction")
                interactions.resolve(debug=debug)
                # self.interactions.resolve(debug=debug, table='temp_interaction')

            ### transfer interactions from temporary table
            self.db.copy_temp_interactions()
            self.set_has_fingerprint(True, commit=commit)

            ### delete temporary table

            if delete_temp_table:
                self.db.execute("DROP TABLE temp_interaction")

        elif debug:
            mrich.warning(f"{self} is already fingerprinted, no new calculation")

    def calculate_prolif_interactions(
        self,
        return_all: bool = False,
        max_retry: int = 5,
        use_mda: bool = False,
        force: bool = False,
        clear_existing: bool = True,
        debug: bool = False,
        resolve: bool = True,
    ) -> "prolif.Fingerprint":
        """Use ProLIF to populate the interactions table"""

        if not self.has_fingerprint or force:

            ### clear old interactions

            if clear_existing:
                self.db.delete_where(
                    table="interaction", key="pose", value=self.id, commit=False
                )
                self.set_has_fingerprint(False, commit=False)

            ### create temporary table

            table = "temp_interaction"

            if "temp_interaction" in self.db.table_names:
                mrich.warning("Deleting existing temp_interaction table")
                self.db.execute("DROP TABLE temp_interaction")

            self.db.create_table_interaction(table="temp_interaction", debug=False)

            if not clear_existing:
                self.db.copy_interactions_to_temp(pose_id=self.id)

            # clear cached InteractionSet
            self._interactions = None

            import prolif as plf
            from tempfile import NamedTemporaryFile
            from .prolif import parse_prolif_interactions
            from MDAnalysis import Universe
            import logging

            mdanalysis_logger = logging.getLogger("MDAnalysis")
            mdanalysis_logger.setLevel(logging.WARNING)

            ## prepare inputs

            # decide if MDA is needed
            unprotonated_sys = self.protein_system
            residue_names = set(r.name for r in unprotonated_sys.residues)
            nonstandard = ["HID", "HIE", "HSE", "HSD", "HSP"]
            if any(r in residue_names for r in nonstandard):
                mrich.debug("Using MDA")
                use_mda = True

            # protonated protein
            for i in range(max_retry):
                try:
                    protonated_sys, protein_file = unprotonated_sys.add_hydrogens(
                        return_file=True
                    )

                    if use_mda:
                        with mrich.loading("Creating MDAnalysis.Universe"):
                            universe = Universe(protein_file.name)
                            protein_mol = plf.Molecule.from_mda(universe)
                    else:
                        with mrich.loading("Creating protein rdkit.Chem.Mol"):
                            rdkit_prot = Chem.MolFromPDBFile(
                                protein_file.name, removeHs=False
                            )
                            protein_mol = plf.Molecule(rdkit_prot)

                    break

                except Exception as e:
                    mrich.warning(
                        f"Could not create satisfactory protein molecule, attempts = {i+1}/{max_retry}"
                    )
                    mrich.warning(e)
                    use_mda = True
                    continue
            else:
                mrich.error(
                    f"Tried {max_retry} times to create protein molecule and failed"
                )
                return None

            # ligand
            ligand_file = NamedTemporaryFile(mode="w+t", suffix=".sdf")
            writer = Chem.SDWriter(ligand_file.name)
            writer.write(self.protonated_mol)
            writer.close()
            ligand_iterable = plf.sdf_supplier(ligand_file.name)

            ## run prolif

            fp = plf.Fingerprint(count=True)
            fp.run_from_iterable(ligand_iterable, protein_mol, progress=False, n_jobs=1)

            ## parse outputs and insert Feature and Interaction records
            parse_prolif_interactions(
                self, fp, protonated_sys, debug=debug, table=table
            )

            if resolve:
                from .iset import InteractionSet

                interactions = InteractionSet.from_pose(self, table="temp_interaction")
                interactions.resolve(debug=debug)

            self.db.copy_temp_interactions()
            self.set_has_fingerprint(True, commit=True)

            ### delete temporary table

            self.db.execute("DROP TABLE temp_interaction")

            ## close files
            protein_file.close()
            ligand_file.close()

            if return_all:
                return fp, ligand_iterable, protein_mol

    def calculate_classic_fingerprint(
        self,
        debug: bool = False,
    ) -> dict:
        """Calculate the pose's interaction fingerprint"""

        if self.path.endswith(".pdb"):

            import molparse as mp

            protein_system = self.protein_system
            if not self.protein_system:
                # mrich.reading(self.path)
                protein_system = mp.parse(self.path, verbosity=False).protein_system

        elif self.path.endswith(".mol") and self.reference:

            # mrich.debug('fingerprint from .mol and reference pose')
            protein_system = self.reference.protein_system

        else:

            mrich.debug("Unsupported: Pose.calculate_fingerprint()")
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
            #   mrich.debug(repr(prot_feature))

            if prot_residue.name != prot_feature.residue_name:
                mrich.warning(f"Feature {repr(prot_feature)}")
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
                    mrich.debug(
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

    def draw2d(
        self,
    ) -> None:
        """Draw a 2D drawing of this pose"""
        from rdkit.Chem import MolFromSmiles

        mol = MolFromSmiles(self.smiles)
        display(mol)

    def render(
        self,
        protein="cartoon",
        ligand="stick",
        protein_color="spectrum",
        interactions: bool = True,
        file: str | None = None,
    ) -> None:
        """Render this pose with the protein using py3Dmol

        :param protein: protein representation, default = 'cartoon'
        :param ligand: ligand representation, default = 'stick'
        :param protein_color: color of protein representation, default = 'spectrum'

        """

        from molparse.py3d import render

        sys = self.complex_system

        def make_view(width="640px", height="480px"):

            view = render(
                sys,
                protein=protein,
                ligand=ligand,
                protein_color=protein_color,
                width=width,
                height=height,
            )

            if interactions:
                COLORS = {
                    "Hydrophobic": "green",
                    "Hydrogen Bond": "blue",
                    "π-stacking": "purple",
                    "π-cation": "pink",
                    "Electrostatic": "red",
                    "Sulfur-Sulfur": "yellow",
                }

                iset = self.interactions

                if not iset:
                    return view

                df = iset.df

                residues = set()

                for i, row in df.iterrows():

                    prot_coord = row["prot_coord"]
                    lig_coord = row["lig_coord"]
                    type = row["type"]
                    color = COLORS.get(type, "black")

                    view.addCylinder(
                        {
                            "start": {
                                "x": prot_coord[0],
                                "y": prot_coord[1],
                                "z": prot_coord[2],
                            },
                            "end": {
                                "x": lig_coord[0],
                                "y": lig_coord[1],
                                "z": lig_coord[2],
                            },
                            # 'radius': radius,
                            "color": color,
                        }
                    )

                    residues.add((row["residue_name"], row["residue_number"]))

                for res_name, res_num in residues:
                    res = sys.residues[f"{res_name} n{res_num}"]
                    view.addModel(res.pdb_block, "pdb")
                    view.setStyle({"model": -1}, {ligand: {}})

            return view

        if file:
            view = make_view(width="100%", height="100%")
            html = view._make_html()

            mrich.writing(file)
            with open(file, "w") as f:
                f.write(html)

        view = make_view()
        return view._repr_html_()

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

    def summary(
        self, metadata: bool = True, tags: bool = True, subsites: bool = True
    ) -> None:
        """Print a summary of this pose

        :param metadata: include metadata (Default value = True)

        """
        if self.alias:
            mrich.header(f"{str(self)}: {self.alias}")
        else:
            mrich.header(f"{str(self)}: {self.inchikey}")
        mrich.var("inchikey", self.inchikey)
        mrich.var("alias", self.alias)
        mrich.var("smiles", self.smiles)
        mrich.var("compound", self.compound)
        mrich.var("path", self.path)
        mrich.var("target", self.target)
        mrich.var("reference", self.reference)
        if tags:
            mrich.var("tags", self.tags)
        if subsites:
            mrich.var("subsites", self.subsites)
        mrich.var("num_heavy_atoms", self.num_heavy_atoms)
        mrich.var("distance_score", self.distance_score)
        mrich.var("energy_score", self.energy_score)
        mrich.var("inspiration_score", self.inspiration_score)
        if inspirations := self.inspirations:
            mrich.var("inspirations", self.inspirations.names)
            mrich.var("num_atoms_added", self.num_atoms_added)
        if metadata:
            mrich.var("metadata", str(self.metadata))

    def showcase(self) -> None:
        """Print and render this pose as if you were using :meth:`.PoseSet.interactive`"""

        self.summary(metadata=False)
        self.grid()
        self.draw()
        from pprint import pprint

        mrich.title("Metadata:")
        pprint(self.metadata)

    def plain_repr(self) -> str:
        """Unformatted detailed string representation"""
        if self.name:
            return f'{self.compound}->{self}: "{self.name}"'
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

    def set_has_fingerprint(self, fp: bool, commit: bool = True) -> None:
        """Update the database to reflect this pose's has_fingerprint property"""
        assert isinstance(fp, bool)
        self._has_fingerprint = fp
        self.db.update(
            table="pose",
            id=self.id,
            key=f"pose_fingerprint",
            value=int(fp),
            commit=commit,
        )

    def posebusters(self, debug: bool = False) -> bool:
        """Run a posebusters ligand check on this pose's molecule"""

        # use syndirella implementation
        from syndirella.slipper import intra_geometry, flatness

        geometries: Dict = intra_geometry.check_geometry(self.mol, threshold_clash=0.4)
        flat_results: Dict = flatness.check_flatness(self.mol)

        if not geometries["results"]["bond_lengths_within_bounds"]:
            if debug:
                mrich.debug(self, "did not pass bond length checks.")
            return False
        if not geometries["results"]["bond_angles_within_bounds"]:
            if debug:
                mrich.debug(self, "did not pass bond angle checks.")
            return False
        if not geometries["results"]["no_internal_clash"]:
            if debug:
                mrich.debug(self, "did not pass internal clash checks.")
            return False
        if not flat_results["results"]["flatness_passes"]:
            if debug:
                mrich.debug(self, "did not pass flatness checks.")
            return False
        return True

    def to_syndirella(self, out_key: "str | Path") -> "DataFrame":
        """Create syndirella inputs. See :meth:`.PoseSet.to_syndirella`"""
        from .pset import PoseSet

        return PoseSet(self.db, [self.id]).to_syndirella(
            out_key=out_key, separate=False
        )

    ### DUNDERS

    def __str__(self) -> str:
        """Unformatted string representation"""
        return f"P{self.id}"

    def __repr__(self) -> str:
        """ANSI Formatted string representation"""
        return f"{mcol.bold}{mcol.underline}{self.plain_repr()}{mcol.unbold}{mcol.ununderline}"

    def __rich__(self) -> str:
        """Formatted string representation"""
        return f"[bold underline]{self.plain_repr()}"

    def __eq__(self, other: "Pose") -> bool:
        """Compare this pose with another instance"""

        if isinstance(other, int):
            return self.id == other

        return self.id == other.id

    def __add__(
        self,
        other: "Pose | PoseSet",
    ) -> "PoseSet":
        """Add a :class:`.PoseSet` to this pose"""
        from .pset import PoseSet

        if isinstance(other, PoseSet):
            return PoseSet(self.db, [self.id] + other.ids, sort=False)
        elif isinstance(other, Pose):
            return PoseSet(self.db, [self.id, other.id], sort=False)
        else:
            raise NotImplementedError


class InvalidMolError(Exception):
    """Exception to be thrown when the molecule could not be parsed"""

    ...
