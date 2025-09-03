import mout
import mcol
import sqlite3
from sqlite3 import Error
from pprint import pprint
import json
import time

from .compound import Compound
from .pose import Pose
from .reaction import Reaction
from .quote import Quote
from .metadata import MetaData
from .target import Target
from .feature import Feature
from .recipe import Recipe, Route
from .tools import inchikey_from_smiles, sanitise_smiles, SanitisationError

from pathlib import Path

import mrich

CHEMICALITE_COMPOUND_PROPERTY_MAP = {
    "num_heavy_atoms": "mol_num_hvyatms",
    "formula": "mol_formula",
    "num_rings": "mol_num_rings",
    "molecular_weight": "mol_amw",
}


class Database:
    """Wrapper to connect to the HIPPO sqlite database.

    .. attention::

            :class:`.Database` objects should not be created directly. Instead use the methods in :class:`.HIPPO` to interact with data in the database. See :doc:`getting_started` and :doc:`insert_elaborations`.

    """

    def __init__(
        self, path: Path, animal: "HIPPO", update_legacy: bool = False
    ) -> None:

        assert isinstance(path, Path)

        mrich.debug("hippo.Database.__init__()")

        self._path = path
        self._connection = None
        self._cursor = None
        self._animal = animal

        mrich.debug(f"Database.path = {self.path}")

        try:
            path = path.resolve(strict=True)

        except FileNotFoundError:
            # create a blank database
            self.connect()
            self.create_blank_db()

        else:
            # connect to existing database
            self.connect()

        if "interaction" not in self.table_names:
            if not update_legacy:
                mrich.error("This is a legacy format database (hippo-db < 0.3.23)")
                mrich.error("Existing fingerprints will not be compatible")
                mrich.error("Re-initialise HIPPO object with update_legacy=True to fix")
                raise LegacyDatabaseError("hippo-db < 0.3.23")
            else:
                mrich.warning("This is a legacy format database (hippo-db < 0.3.23)")
                mrich.warning("Clearing legacy fingerprints...")
                self.create_table_interaction()
                self.delete_interactions()

        if "subsite" not in self.table_names or "subsite_tag" not in self.table_names:
            mrich.warning("This is a legacy format database (hippo-db < 0.3.24)")
            self.create_table_subsite()
            self.create_table_subsite_tag()

        if "scaffold" not in self.table_names:
            if not update_legacy:
                mrich.error("This is a legacy format database (hippo-db < 0.3.25)")
                mrich.error("Existing base-elab relationships will not be compatible")
                mrich.error("Re-initialise HIPPO object with update_legacy=True to fix")
                raise LegacyDatabaseError("hippo-db < 0.3.25")
            else:
                mrich.warning("This is a legacy format database (hippo-db < 0.3.25)")
                mrich.warning("Migrating compound_base values to scaffold table...")
                self.create_table_scaffold()
                self.migrate_legacy_bases()

        if "route" not in self.table_names:
            self.create_table_route()
            self.create_table_component()

        elif "component_amount" not in self.column_names("component"):
            if not update_legacy:
                mrich.error("This is a legacy format database (hippo-db < 0.3.29)")
                mrich.error("Re-initialise HIPPO object with update_legacy=True to fix")
                raise LegacyDatabaseError("hippo-db < 0.3.29")
            else:
                mrich.warning("This is a legacy format database (hippo-db < 0.3.29)")
                mrich.warning("Updating legacy routes...")

            self.update_legacy_routes()

        if "reaction_metadata" not in self.column_names("reaction"):
            if not update_legacy:
                mrich.error("This is a legacy format database (hippo-db < 0.3.32)")
                mrich.error("Re-initialise HIPPO object with update_legacy=True to fix")
                raise LegacyDatabaseError("hippo-db < 0.3.32")
            else:
                mrich.warning("This is a legacy format database (hippo-db < 0.3.32)")
                mrich.warning("Updating legacy reaction table...")

            self.update_legacy_reaction_metadata()

    @classmethod
    def copy_from(
        cls,
        source: Path,
        destination: Path,
        animal: "HIPPO",
        update_legacy: bool = False,
        overwrite_existing: bool = False,
        pages: int = 10000,
    ) -> None:
        """Create a :class:`.Database` from an existing one"""

        source = Path(source)

        assert source.exists()

        if destination.exists():
            if overwrite_existing:
                mrich.warning(f"Overwriting {destination}")
            else:
                mrich.error(f"Will not overwrite {destination}")
                mrich.error(f"Set overwrite_existing=True to override")
                raise Exception("Set overwrite_existing=True to override")

        mrich.print(f"Copying {source} --> {destination}")

        def progress(status, remaining, total):
            mrich.debug(f"Copied {total-remaining} of {total} pages...")

        src = sqlite3.connect(source)
        dst = sqlite3.connect(destination)
        with dst:
            src.backup(dst, pages=pages, progress=progress)
        dst.close()
        src.close()

        self = cls.__new__(cls)

        self.__init__(path=destination, animal=animal, update_legacy=update_legacy)

        return self

    ### PROPERTIES

    @property
    def path(self) -> Path:
        """Returns the path to the database file"""
        return self._path

    @property
    def connection(self) -> "sqlite3.connection":
        """Returns a ``sqlite3.connection`` to the database"""
        if not self._connection:
            self.connect()
        return self._connection

    @property
    def cursor(self) -> "sqlite3.cursor":
        """Returns a ``sqlite3.cursor``"""
        if not self._cursor:
            self._cursor = self.connection.cursor()
        return self._cursor

    @property
    def total_changes(self) -> int:
        """Return the total number of database rows that have been modified, inserted, or deleted since the database connection was opened."""
        return self.connection.total_changes

    @property
    def table_names(self) -> list[str]:
        """List of all the table names in the database"""
        results = self.execute(
            "SELECT name FROM sqlite_master WHERE type='table';"
        ).fetchall()
        return [n for n, in results]

    ### PUBLIC METHODS / API CALLS

    def close(self) -> None:
        """Close the connection"""
        mrich.debug("hippo.Database.close()")
        if self.connection:
            self.connection.close()
        mout.success(f"Closed connection to {self.path}")

    def backup(
        self, destination: Path | str | None = None, pages: int = 10_000
    ) -> None:
        """Create a backup of the database"""

        from .tools import dt_hash

        if not destination:
            destination = str(self.path.resolve()).replace(
                ".sqlite", f"_{dt_hash()}.sqlite"
            )

        destination = Path(destination)

        with mrich.spinner(f"Backing up {self.__rich__()}[reset]..."):

            mrich.writing(destination)

            def progress(status, remaining, total):
                mrich.debug(f"Copied {total-remaining} of {total} pages...")

            src = self.connection
            dst = sqlite3.connect(destination)
            with dst:
                src.backup(dst, pages=pages, progress=progress)

            dst.close()

    ### GENERAL SQL

    def connect(self) -> None:
        """Connect to the database"""
        mrich.debug("hippo.Database.connect()")

        conn = None

        try:
            conn = sqlite3.connect(self.path)

            mrich.debug(f"{sqlite3.version=}")

            conn.enable_load_extension(True)
            conn.load_extension("chemicalite")
            conn.enable_load_extension(False)

            mrich.success("Database connected @", f"[file]{self.path}")

        except sqlite3.OperationalError as e:

            if "cannot open shared object file" in str(e):
                mrich.error("chemicalite package not installed correctly")
            else:
                mrich.error(e)
            raise

        except Error as e:
            mrich.error(e)
            raise

        self._connection = conn
        self._cursor = conn.cursor()

    def execute(
        self, sql, payload=None, *, retry: float | None = 1, debug: bool = False
    ) -> None:
        """Execute arbitrary SQL

        :param sql: SQL query
        :param retry: If truthy, keep trying to execute every `retry` seconds if the Database is locked
        :param payload: Payload for insertion, etc. (Default value = None)

        """

        if debug:
            mrich.debug(sql)

        try:
            if payload:
                return self.cursor.execute(sql, payload)
            else:
                return self.cursor.execute(sql)
        except sqlite3.OperationalError as e:
            if "database is locked" in str(e) and retry:
                with mrich.clock(
                    f"SQLite Database is locked, waiting {retry} second(s)..."
                ):
                    time.sleep(retry)
                mrich.print("[debug]SQLite Database was locked, retrying...")
                return self.execute(sql=sql, payload=payload, retry=retry)
            if "syntax error" in str(e):
                mrich.error(sql)
                mrich.error(payload)
                raise
            else:
                raise
        except Exception as e:
            raise

    def executemany(self, sql, payload, *, retry: float | None = 1) -> None:
        """Execute arbitrary SQL

        :param sql: SQL query
        :param retry: If truthy, keep trying to executemany every `retry` seconds if the Database is locked
        :param payload: Payload for insertion, etc. (Default value = None)

        """

        if "RETURNING" in sql:
            from .apsw import executemany

            return executemany(self.path, sql, payload)

        try:
            return self.cursor.executemany(sql, payload)
        except sqlite3.OperationalError as e:
            if "database is locked" in str(e) and retry:
                with mrich.clock(
                    f"SQLite Database is locked, waiting {retry} second(s)..."
                ):
                    time.sleep(retry)
                mrich.print("[debug]SQLite Database was locked, retrying...")
                return self.executemany(sql=sql, payload=payload, retry=retry)
            else:
                raise
        except Exception as e:
            mrich.print(sql)
            mrich.print(payload[0])
            raise

    def commit(self, *, retry: float | None = 1) -> None:
        """Commit changes to the database

        :param retry: If truthy, keep trying to execute every `retry` seconds if the Database is locked
        """
        try:
            self.connection.commit()
        except sqlite3.OperationalError as e:
            if "database is locked" in str(e) and retry:
                with mrich.clock(
                    f"SQLite Database is locked, waiting {retry} second(s)..."
                ):
                    time.sleep(retry)
                mrich.print("[debug]SQLite Database was locked, retrying...")
                return self.commit()
            else:
                raise

    ### CREATE TABLES

    def create_blank_db(self) -> None:
        """Create a blank database"""

        with mrich.loading("Creating blank database..."):
            self.create_table_compound()
            self.create_table_inspiration()
            self.create_table_reaction()
            self.create_table_reactant()
            self.create_table_pose()
            self.create_table_tag()
            self.create_table_quote()
            self.create_table_target()
            self.create_table_pattern_bfp()
            self.create_table_feature()
            self.create_table_route()
            self.create_table_component()
            self.create_table_interaction()
            self.create_table_subsite()
            self.create_table_subsite_tag()
            self.create_table_scaffold()
            self.commit()

    def create_table_compound(self) -> None:
        """Create the compound table"""
        mrich.debug("HIPPO.Database.create_table_compound()")

        sql = """CREATE TABLE compound(
            compound_id INTEGER PRIMARY KEY,
            compound_inchikey TEXT,
            compound_alias TEXT,
            compound_smiles TEXT,
            compound_base INTEGER,
            compound_mol MOL,
            compound_pattern_bfp bits(2048),
            compound_morgan_bfp bits(2048),
            compound_metadata TEXT,
            FOREIGN KEY (compound_base) REFERENCES compound(compound_id),
            CONSTRAINT UC_compound_inchikey UNIQUE (compound_inchikey)
            CONSTRAINT UC_compound_alias UNIQUE (compound_alias)
            CONSTRAINT UC_compound_smiles UNIQUE (compound_smiles)
        );
        """
        self.execute(sql)

    def create_table_inspiration(self) -> None:
        """Create the inspiration table"""
        mrich.debug("HIPPO.Database.create_table_inspiration()")

        sql = """CREATE TABLE inspiration(
            inspiration_original INTEGER,
            inspiration_derivative INTEGER,
            FOREIGN KEY (inspiration_original) REFERENCES pose(pose_id),
            FOREIGN KEY (inspiration_derivative) REFERENCES pose(pose_id),
            CONSTRAINT UC_inspiration UNIQUE (inspiration_original, inspiration_derivative)
        );
        """

        self.execute(sql)

    def create_table_scaffold(self) -> None:
        """Create the scaffold table"""
        mrich.debug("HIPPO.Database.create_table_scaffold()")

        sql = """CREATE TABLE scaffold(
            scaffold_base INTEGER,
            scaffold_superstructure INTEGER,
            FOREIGN KEY (scaffold_base) REFERENCES pose(pose_id),
            FOREIGN KEY (scaffold_superstructure) REFERENCES pose(pose_id),
            CONSTRAINT UC_scaffold UNIQUE (scaffold_base, scaffold_superstructure)
        );
        """

        self.execute(sql)

    def create_table_reaction(self) -> None:
        """Create the reaction table"""
        mrich.debug("HIPPO.Database.create_table_reaction()")

        sql = """CREATE TABLE reaction(
            reaction_id INTEGER PRIMARY KEY,
            reaction_type TEXT,
            reaction_product INTEGER,
            reaction_product_yield REAL,
            reaction_metadata TEXT,
            FOREIGN KEY (reaction_product) REFERENCES compound(compound_id)
        );
        """

        self.execute(sql)

    def create_table_reactant(self) -> None:
        """Create the reactant table"""
        mrich.debug("HIPPO.Database.create_table_reactant()")

        sql = """CREATE TABLE reactant(
            reactant_amount REAL,
            reactant_reaction INTEGER,
            reactant_compound INTEGER,
            FOREIGN KEY (reactant_reaction) REFERENCES reaction(reaction_id),
            FOREIGN KEY (reactant_compound) REFERENCES compound(compound_id),
            CONSTRAINT UC_reactant UNIQUE (reactant_reaction, reactant_compound)
        );
        """

        self.execute(sql)

    def create_table_pose(self) -> None:
        """Create the pose table"""
        mrich.debug("HIPPO.Database.create_table_pose()")

        sql = """CREATE TABLE pose(
            pose_id INTEGER PRIMARY KEY,
            pose_inchikey TEXT,
            pose_alias TEXT,
            pose_smiles TEXT,
            pose_reference INTEGER,
            pose_path TEXT,
            pose_compound INTEGER,
            pose_target INTEGER,
            pose_mol BLOB,
            pose_fingerprint BLOB,
            pose_energy_score REAL,
            pose_distance_score REAL,
            pose_metadata TEXT,
            FOREIGN KEY (pose_compound) REFERENCES compound(compound_id),
            CONSTRAINT UC_pose_alias UNIQUE (pose_alias)
            CONSTRAINT UC_pose_path UNIQUE (pose_path)
        );
        """

        ### snippet to convert python metadata dictionary with JSON
        # json.dumps(variables).encode('utf-8')
        # json.loads(s.decode('utf-8'))

        self.execute(sql)

    def create_table_tag(self) -> None:
        """Create the tag table"""
        mrich.debug("HIPPO.Database.create_table_tag()")

        sql = """CREATE TABLE tag(
            tag_name TEXT,
            tag_compound INTEGER,
            tag_pose INTEGER,
            FOREIGN KEY (tag_compound) REFERENCES compound(compound_id),
            FOREIGN KEY (tag_pose) REFERENCES pose(pose_id),
            CONSTRAINT UC_tag_compound UNIQUE (tag_name, tag_compound)
            CONSTRAINT UC_tag_pose UNIQUE (tag_name, tag_pose)
        );
        """

        self.execute(sql)

    def create_table_quote(self) -> None:
        """Create the quote table"""
        mrich.debug("HIPPO.Database.create_table_quote()")

        sql = """CREATE TABLE quote(
            quote_id INTEGER PRIMARY KEY,
            quote_smiles TEXT,
            quote_amount REAL,
            quote_supplier TEXT,
            quote_catalogue TEXT,
            quote_entry TEXT,
            quote_lead_time INTEGER,
            quote_price REAL,
            quote_currency TEXT,
            quote_purity REAL,
            quote_date TEXT,
            quote_compound INTEGER,
            FOREIGN KEY (quote_compound) REFERENCES compound(compound_id),
            CONSTRAINT UC_quote UNIQUE (quote_amount, quote_supplier, quote_catalogue, quote_entry)
        );
        """

        self.execute(sql)

    def create_table_target(self) -> None:
        """Create the target table"""
        mrich.debug("HIPPO.Database.create_table_target()")
        sql = """CREATE TABLE target(
            target_id INTEGER PRIMARY KEY,
            target_name TEXT,
            target_metadata TEXT,
            CONSTRAINT UC_target UNIQUE (target_name)
        );
        """

        self.execute(sql)

    def create_table_feature(self) -> None:
        """Create the feature table"""
        mrich.debug("HIPPO.Database.create_table_feature()")
        sql = """CREATE TABLE feature(
            feature_id INTEGER PRIMARY KEY,
            feature_family TEXT,
            feature_target INTEGER,
            feature_chain_name TEXT,
            feature_residue_name TEXT,
            feature_residue_number INTEGER,
            feature_atom_names TEXT,
            CONSTRAINT UC_feature UNIQUE (feature_family, feature_target, feature_chain_name, feature_residue_number, feature_residue_name, feature_atom_names)
        );
        """

        self.execute(sql)

    def create_table_route(self) -> None:
        """Create the route table"""
        mrich.debug("HIPPO.Database.create_table_route()")
        sql = """CREATE TABLE route(
            route_id INTEGER PRIMARY KEY,
            route_product INTEGER
        );
        """

        self.execute(sql)

    def create_table_component(self) -> None:
        """Create the component table"""
        mrich.debug("HIPPO.Database.create_table_component()")
        sql = """CREATE TABLE component(
            component_id INTEGER PRIMARY KEY,
            component_route INTEGER,
            component_type INTEGER,
            component_ref INTEGER,
            component_amount REAL,
            CONSTRAINT UC_component UNIQUE (component_route, component_ref, component_type)
        );
        """

        self.execute(sql)

    def create_table_pattern_bfp(self) -> None:
        """Create the pattern_bfp table"""
        mrich.debug("HIPPO.Database.create_table_pattern_bfp()")

        sql = "CREATE VIRTUAL TABLE compound_pattern_bfp USING rdtree(compound_id, fp bits(2048))"

        self.execute(sql)

    def create_table_interaction(
        self, table: str = "interaction", debug: bool = True
    ) -> None:
        """Create an interaction table"""

        if debug:
            mrich.debug(f"HIPPO.Database.create_table_interaction({table=})")
        sql = f"""CREATE TABLE {table}(
            interaction_id INTEGER PRIMARY KEY,
            interaction_feature INTEGER NOT NULL,
            interaction_pose INTEGER NOT NULL,
            interaction_type TEXT NOT NULL,
            interaction_family TEXT NOT NULL,
            interaction_atom_ids TEXT NOT NULL,
            interaction_prot_coord TEXT NOT NULL,
            interaction_lig_coord TEXT NOT NULL,
            interaction_distance REAL NOT NULL,
            interaction_angle REAL,
            interaction_energy REAL,
            CONSTRAINT UC_interaction UNIQUE (interaction_feature, interaction_pose, interaction_type, interaction_family, interaction_atom_ids)
        );
        """

        self.execute(sql)

    def create_table_subsite(self) -> None:
        """Create the subsite table"""

        mrich.debug("HIPPO.Database.create_table_subsite()")
        sql = """CREATE TABLE subsite(
            subsite_id INTEGER PRIMARY KEY,
            subsite_target INTEGER NOT NULL,
            subsite_name TEXT NOT NULL,
            subsite_metadata TEXT,
            CONSTRAINT UC_subsite UNIQUE (subsite_target, subsite_name)
        );
        """

        self.execute(sql)

    def create_table_subsite_tag(self) -> None:
        """Create the subsite_tag table"""

        mrich.debug("HIPPO.Database.create_table_subsite_tag()")
        sql = """CREATE TABLE subsite_tag(
            subsite_tag_id INTEGER PRIMARY KEY,
            subsite_tag_ref INTEGER NOT NULL,
            subsite_tag_pose INTEGER NOT NULL,
            subsite_tag_metadata TEXT,
            CONSTRAINT UC_subsite_tag UNIQUE (subsite_tag_ref, subsite_tag_pose)
        );
        """

        self.execute(sql)

    ### INSERTION

    def insert_compound(
        self,
        *,
        smiles: str,
        alias: str | None = None,
        tags: None | list[str] = None,
        warn_duplicate: bool = True,
        commit: bool = True,
        metadata: None | dict = None,
        inchikey: str = None,
    ) -> int:
        """Insert an entry into the compound table

        :param smiles: SMILES string
        :param alias: optional alias for the compound (Default value = None)
        :param tags: list of string tags, (Default value = None)
        :param warn_duplicate: print a warning if the compound already exists (Default value = True)
        :param commit: commit the changes to the database (Default value = True)
        :param metadata: dictionary of metadata (Default value = None)
        :param inchikey: provide an InChI-key, otherwise it's calculated from the SMILES, (Default value = None)
        :returns: compound ID

        """

        # generate the inchikey name
        inchikey = inchikey or inchikey_from_smiles(smiles)

        sql = """
        INSERT INTO compound(compound_inchikey, compound_smiles, compound_mol, compound_pattern_bfp, compound_morgan_bfp, compound_alias)
        VALUES(?1, ?2, mol_from_smiles(?2), mol_pattern_bfp(mol_from_smiles(?2), 2048), mol_morgan_bfp(mol_from_smiles(?2), 2, 2048), ?3)
        """

        try:
            self.execute(sql, (inchikey, smiles, alias))

        except sqlite3.IntegrityError as e:
            if "UNIQUE constraint failed: compound.compound_inchikey" in str(e):
                if warn_duplicate:
                    mrich.warning(
                        f'Skipping compound with existing inchikey "{inchikey}"'
                    )
            elif "UNIQUE constraint failed: compound.compound_smiles" in str(e):
                if warn_duplicate:
                    mrich.warning(f'Skipping compound with existing smiles "{smiles}"')
            elif "UNIQUE constraint failed: compound.compound_pattern_bfp" in str(e):
                if warn_duplicate:
                    mrich.warning(
                        f'Skipping compound with existing pattern binary fingerprint "{smiles}"'
                    )
            elif "UNIQUE constraint failed: compound.compound_morgan_bfp" in str(e):
                if warn_duplicate:
                    mrich.warning(
                        f'Skipping compound with existing morgan binary fingerprint "{smiles}"'
                    )
            else:
                mrich.error(e)

            return None

        except Exception as e:
            mrich.error(e)

        compound_id = self.cursor.lastrowid
        if commit:
            self.commit()

        ### register the tags

        if tags:
            for tag in tags:
                self.insert_tag(name=tag, compound=compound_id, commit=commit)

        ### register the binary fingerprints

        result = self.insert_compound_pattern_bfp(compound_id, commit=commit)

        if not result:
            mrich.error("Could not insert compound pattern bfp")

        if metadata:
            self.insert_metadata(
                table="compound", id=compound_id, payload=metadata, commit=commit
            )

        return compound_id

    def insert_compound_pattern_bfp(self, compound_id: int, commit: bool = True) -> int:
        """Insert a compound_pattern_bfp

        :param compound_id: ID of the associated compound
        :param commit: commit the changes to the database (Default value = True)
        :returns: binary fingerprint ID

        """

        sql = """
        INSERT INTO compound_pattern_bfp(compound_id, fp)
        VALUES(?1, ?2)
        """

        (bfp,) = self.select_where(
            "compound_pattern_bfp", "compound", "id", compound_id
        )

        try:
            self.execute(sql, (compound_id, bfp))

        except Exception as e:
            mrich.error(e)

        bfp_id = self.cursor.lastrowid
        if commit:
            self.commit()

        return bfp_id

    def insert_pose(
        self,
        *,
        compound: Compound | int,
        target: Target | int | str,
        path: str,
        inchikey: str | None = None,
        alias: str | None = None,
        reference: int | Pose | None = None,
        tags: None | list = None,
        energy_score: float | None = None,
        distance_score: float | None = None,
        metadata: None | dict = None,
        commit: bool = True,
        warn_duplicate: bool = True,
        resolve_path: bool = True,
    ) -> int:
        """Insert an entry into the pose table

        :param compound: associated :class:`.Compound` object or ID
        :param target: protein :class:`.Target` name or ID
        :param path: path to the molecular structure (.pdb/.mol)
        :param inchikey: provide an InChI-key if available, (Default value = None)
        :param alias: optional alias for the compound (Default value = None)
        :param reference: reference :class:`.Pose` object or ID to use for the protein conformation (Default value = None)
        :param tags: list of string tags, (Default value = None)
        :param energy_score: optional score of the ligand's binding energy (Default value = None)
        :param distance_score: optional score of the ligand's binding position (Default value = None)
        :param metadata: dictionary of metadata (Default value = None)
        :param commit: commit the changes to the database (Default value = True)
        :param warn_duplicate: print a warning if the pose already exists (Default value = True)
        :param resolve_path: try resolving the path (Default value = True)
        :returns: the pose ID

        """

        if isinstance(compound, Compound):
            compound = compound.id

        if isinstance(reference, Pose):
            reference = reference.id

        if isinstance(target, Target):
            target = target.id

        if isinstance(target, str):
            target = self.get_target_id(name=target)
            if not target:
                raise ValueError(f"No such {target=}")
        target_name = self.get_target_name(id=target)

        if resolve_path:
            try:
                path = Path(path)
                path = path.resolve(strict=True)
                path = str(path)

            except FileNotFoundError as e:
                mrich.error(f"Path cannot be resolved: {mcol.file}{path}")
                raise

        sql = """
        INSERT INTO pose(pose_inchikey, pose_alias, pose_smiles, pose_compound, pose_target, pose_path, pose_reference, pose_energy_score, pose_distance_score)
        VALUES(?1, ?2, ?3, ?4, ?5, ?6, ?7, ?8, ?9)
        """

        try:
            self.execute(
                sql,
                (
                    inchikey,
                    alias,
                    None,
                    compound,
                    target,
                    path,
                    reference,
                    energy_score,
                    distance_score,
                ),
            )

        except sqlite3.IntegrityError as e:
            if "UNIQUE constraint failed: pose.pose_path" in str(e):
                if warn_duplicate:
                    mrich.warning(f'Could not insert pose with duplicate path "{path}"')
            elif "UNIQUE constraint failed: pose.pose_alias" in str(e):
                if warn_duplicate:
                    mrich.warning(
                        f'Could not insert pose with duplicate alias "{alias}"'
                    )
            else:
                mrich.error(e)
            return None

        except Exception as e:
            mrich.error(e)
            raise

        pose_id = self.cursor.lastrowid
        if commit:
            self.commit()

        if tags:
            for tag in tags:
                self.insert_tag(name=tag, pose=pose_id, commit=commit)

        if metadata:
            self.insert_metadata(
                table="pose", id=pose_id, payload=metadata, commit=commit
            )

        return pose_id

    def insert_tag(
        self,
        *,
        name: str,
        compound: int = None,
        pose: int = None,
        commit: bool = True,
    ) -> int:
        """Insert an entry into the tag table.

        .. attention::
                Exactly one of compound or pose arguments must have a value

        :param name: name of the tag
        :param compound: associated :class:`.Compound` ID
        :param pose: associated :class:`.Pose` ID
        :param commit: commit the changes to the database (Default value = True)
        :returns: the tag ID

        """

        assert bool(compound) ^ bool(
            pose
        ), "Exactly one of compound or pose arguments must have a value"

        sql = """
        INSERT INTO tag(tag_name, tag_compound, tag_pose)
        VALUES(?1, ?2, ?3)
        """

        try:
            self.execute(sql, (name, compound, pose))

        except sqlite3.IntegrityError as e:
            return None

        except Exception as e:
            mrich.error(e)

        tag_id = self.cursor.lastrowid
        if commit:
            self.commit()
        return tag_id

    def insert_inspiration(
        self,
        *,
        original: Pose | int,
        derivative: Pose | int,
        warn_duplicate: bool = True,
        commit: bool = True,
    ) -> int:
        """Insert an entry into the inspiration table

        :param original: :class:`.Pose` object or ID of the original hit
        :param derivative: :class:`.Pose` object or ID of the derivative hit
        :param warn_duplicate: print a warning if the pose already exists (Default value = True)
        :param commit: commit the changes to the database (Default value = True)
        :returns: the inspiration ID

        """

        if isinstance(original, Pose):
            original = original.id
        if isinstance(derivative, Pose):
            derivative = derivative.id

        assert isinstance(
            original, int
        ), "Must pass an integer ID or Pose object (original)"
        assert isinstance(
            derivative, int
        ), "Must pass an integer ID or Pose object (derivative)"

        sql = """
        INSERT INTO inspiration(inspiration_original, inspiration_derivative)
        VALUES(?1, ?2)
        """

        try:
            self.execute(sql, (original, derivative))

        except sqlite3.IntegrityError as e:
            if warn_duplicate:
                mrich.warning(
                    f"Skipping existing inspiration: {original=} {derivative=}"
                )
            return None

        except Exception as e:
            mrich.error(e)

        inspiration_id = self.cursor.lastrowid

        if commit:
            self.commit()
        return inspiration_id

    def insert_scaffold(
        self,
        *,
        base: Compound | int,
        superstructure: Compound | int,
        warn_duplicate: bool = True,
        commit: bool = True,
    ) -> int:
        """Insert an entry into the scaffold table

        :param base: :class:`.Compound` object or ID of the base hit
        :param superstructure: :class:`.Compound` object or ID of the superstructure hit
        :param warn_duplicate: print a warning if the pose already exists (Default value = True)
        :param commit: commit the changes to the database (Default value = True)
        :returns: the scaffold row ID

        """

        if isinstance(base, Compound):
            base = base.id
        if isinstance(superstructure, Compound):
            superstructure = superstructure.id

        assert isinstance(
            base, int
        ), f"Must pass an integer ID or Compound object (base) {base=} {type(base)}"
        assert isinstance(
            superstructure, int
        ), f"Must pass an integer ID or Compound object (superstructure) {superstructure=} {type(superstructure)}"

        if base == superstructure:
            # mrich.warning(f"Skipped self-referential scaffold assignment (C{base})")
            return None

        sql = """
        INSERT INTO scaffold(scaffold_base, scaffold_superstructure)
        VALUES(?1, ?2)
        """

        try:
            self.execute(sql, (base, superstructure))

        except sqlite3.IntegrityError as e:
            if warn_duplicate:
                mrich.warning(f"Skipping existing scaffold: {base=} {superstructure=}")
            return None

        except Exception as e:
            mrich.error(e)

        scaffold_id = self.cursor.lastrowid

        if commit:
            self.commit()
        return scaffold_id

    def insert_reaction(
        self,
        *,
        type: str,
        product: Compound | int,
        product_yield: float = 1.0,
        commit: bool = True,
    ) -> int:
        """Insert an entry into the reaction table

        :param type: string to indicate the reaction type
        :param product: :class:`.Compound` object or ID of the reaction product
        :param product_yield: yield fraction of the reaction product (Default value = 1.0)
        :param commit: commit the changes to the database (Default value = True)
        :returns: the reaction ID

        """

        if isinstance(product, Compound):
            product = product.id

        # assert isinstance(product, Compound), f'incompatible {product=}'
        assert isinstance(type, str), f"incompatible {type=}"

        sql = """
        INSERT INTO reaction(reaction_type, reaction_product, reaction_product_yield)
        VALUES(?1, ?2, ?3)
        """

        try:
            self.execute(sql, (type, product, product_yield))

        except Exception as e:
            mrich.error(e)

        reaction_id = self.cursor.lastrowid
        if commit:
            self.commit()
        return reaction_id

    def insert_reactant(
        self,
        *,
        compound: Compound | int,
        reaction: Reaction | int,
        amount: float = 1.0,
        commit: bool = True,
    ) -> int:
        """Insert an entry into the reactant table

        :param compound: :class:`.Compound` object or ID of the reactant
        :param reaction: :class:`.Reaction` object or ID of the reaction
        :param amount: amount (in ``mg``) needed for each unit of product (Default value = 1.0)
        :param commit: commit the changes to the database (Default value = True)
        :returns: the reactant ID

        """

        if isinstance(reaction, int):
            reaction = self.get_reaction(id=reaction)

        if isinstance(compound, int):
            compound = self.get_compound(id=compound)

        assert isinstance(compound, Compound), f"incompatible {compound=}"
        assert isinstance(reaction, Reaction), f"incompatible {reaction=}"

        sql = """
        INSERT INTO reactant(reactant_amount, reactant_reaction, reactant_compound)
        VALUES(?1, ?2, ?3)
        """

        try:
            self.execute(sql, (amount, reaction.id, compound.id))

        except sqlite3.IntegrityError as e:
            mrich.warning(f"Skipping existing reactant: {reaction=} {compound=}")

        except Exception as e:
            mrich.error(e)

        reactant_id = self.cursor.lastrowid

        if commit:
            self.commit()

        return reactant_id

    def insert_quote(
        self,
        *,
        compound: Compound | int,
        supplier: str,
        catalogue: str | None = None,
        entry: str | None = None,
        amount: float,
        price: float,
        currency: str | None = None,
        purity: float | None = None,
        lead_time: float,
        smiles: str | None = None,
        date: str | None = None,
        commit: bool = True,
    ) -> int | None:
        """Insert an entry into the quote table

        :param compound: associated :class:`.Compound` object or ID
        :param supplier: name of the supplier
        :param catalogue: optional catalogue name
        :param entry: name of the catalogue entry
        :param amount: amount in `mg`
        :param price: price of the compound
        :param currency: currency string ``['GBP', 'EUR', 'USD', None]``
        :param purity: compound purity fraction
        :param lead_time: lead time in days
        :param smiles: quoted SMILES string (Default value = None)
        :param commit: commit the changes to the database (Default value = True)
        :returns: the quote ID

        """

        if not isinstance(compound, int):
            assert isinstance(compound, Compound), f"incompatible {compound=}"
            compound = compound.id

        assert currency in ["GBP", "EUR", "USD", None], f"incompatible {currency=}"

        assert supplier in [
            "MCule",
            "Enamine",
            "Stock",
            "Molport",
        ], f"incompatible {supplier=}"

        smiles = smiles or ""

        payload = [
            smiles,
            amount,
            supplier,
            catalogue,
            entry,
            lead_time,
            price,
            currency,
            purity,
            compound,
        ]

        if date:
            date_str = "?11"
            payload.append(date)
        else:
            date_str = "date()"

        sql = f"""
        INSERT or REPLACE INTO quote(
            quote_smiles,
            quote_amount,
            quote_supplier,
            quote_catalogue,
            quote_entry,
            quote_lead_time,
            quote_price,
            quote_currency,
            quote_purity,
            quote_compound,
            quote_date
        )
        VALUES(?1, ?2, ?3, ?4, ?5, ?6, ?7, ?8, ?9, ?10, {date_str});
        """

        try:
            self.execute(
                sql,
                tuple(payload),
            )

        except sqlite3.InterfaceError as e:
            mrich.error(e)
            mrich.debug(payload)
            raise

        except Exception as e:
            mrich.error(e)
            return None

        quote_id = self.cursor.lastrowid
        if commit:
            self.commit()
        return quote_id

    def insert_target(
        self,
        *,
        name: str,
        warn_duplicate: bool = True,
    ) -> int:
        """Insert an entry into the target table

        :param name: name of the protein target
        :returns: the target ID

        """

        sql = """
        INSERT INTO target(target_name)
        VALUES(?1)
        """

        try:
            self.execute(sql, (name,))

        except sqlite3.IntegrityError as e:
            if warn_duplicate:
                mrich.warning(f"Skipping existing target with {name=}")
            return None

        except Exception as e:
            mrich.error(e)

        target_id = self.cursor.lastrowid
        self.commit()
        return target_id

    def insert_feature(
        self,
        *,
        family: str,
        target: int,
        chain_name: str,
        residue_name: str,
        residue_number: int,
        atom_names: list[str],
        warn_duplicate: bool = False,
        commit: bool = True,
    ) -> int:
        """Insert an entry into the feature table

        :param family: feature type string
        :param target: associated :class:`.Target` ID
        :param chain_name: single character name of the chain
        :param residue_name: 3-4 character string name of the residue
        :param residue_number: integer residue number
        :param atom_names: list of atom names
        :param commit: commit the changes to the database (Default value = True)
        :returns: feature ID

        """

        assert len(chain_name) == 1
        assert len(residue_name) <= 4
        for a in atom_names:
            assert len(a) <= 4

        if isinstance(target, str):
            target = self.get_target_id(name=target)
        assert isinstance(target, int)

        from .prolif import FEATURE_FAMILIES

        if family:
            assert family in FEATURE_FAMILIES, f"Unsupported {family=}"
        else:
            family = "Unknown"

        sql = """
        INSERT INTO feature(feature_family, feature_target, feature_chain_name, feature_residue_name, feature_residue_number, feature_atom_names)
        VALUES(?1, ?2, ?3, ?4, ?5, ?6)
        """

        atom_names = " ".join(sorted(atom_names))

        try:
            self.execute(
                sql,
                (family, target, chain_name, residue_name, residue_number, atom_names),
            )

        except sqlite3.IntegrityError as e:

            if warn_duplicate:
                mrich.warning(str(e))
                mrich.var("family", family)
                mrich.var("target", target)
                mrich.var("chain_name", chain_name)
                mrich.var("residue_name", residue_name)
                mrich.var("residue_number", residue_number)
                mrich.var("atom_names", atom_names)

            return None

        except Exception as e:
            mrich.error(e)

        feature_id = self.cursor.lastrowid
        if commit:
            self.commit()
        return feature_id

    def insert_metadata(
        self,
        *,
        table: str,
        id: int,
        payload: dict,
        commit: bool = True,
    ) -> None:
        """Insert metadata into an an existing entry in the compound or pose tables

        :param table: table for insertions ``['pose', 'compound', 'subsite', 'subsite_tag']``
        :param id: associated entry ID
        :param payload: metadata dictionary
        :param commit: commit the changes to the database (Default value = True)

        """

        payload = json.dumps(payload)

        self.update(
            table=table, id=id, key=f"{table}_metadata", value=payload, commit=commit
        )

    def insert_route(
        self,
        *,
        product_id: int,
        commit: bool = True,
    ) -> int:
        """Insert an entry into the route table

        :param product_id: :class:`.Compound` ID of the product
        :param commit: commit the changes to the database (Default value = True)
        :returns: route ID

        """

        sql = """
        INSERT INTO route(route_product)
        VALUES(?1)
        """

        product_id = int(product_id)

        try:
            self.execute(sql, (product_id,))

        except Exception as e:
            mrich.error(e)

        route_id = self.cursor.lastrowid

        if commit:
            self.commit()

        return route_id

    def insert_component(
        self,
        *,
        route: int,
        ref: int,
        component_type: int,
        amount: float = 1.0,
        commit: bool = True,
    ) -> int:
        """
        ================ ========== ==============
        component_type   table      ref
        ================ ========== ==============
                      1   reaction   reaction
                      2   compound   reactant
                      3   compound   intermediate
        ================ ========== ==============

        :param route: associated :class:`.Route` ID
        :param ref: ID of the :class:`.Reaction` or :class:`.Compound`
        :param component_type: integer specifying the type of the component
        :param commit: commit the changes to the database (Default value = True)
        :returns: the component ID

        """

        sql = """
        INSERT INTO component(component_route, component_type, component_ref, component_amount)
        VALUES(:component_route, :component_type, :component_ref, :component_amount)
        """

        route = int(route)
        ref = int(ref)
        component_type = int(component_type)

        if component_type == 1:
            component_amount = None
        else:
            component_amount = float(amount)
            assert component_amount > 0

        try:
            self.execute(
                sql,
                dict(
                    component_route=route,
                    component_type=component_type,
                    component_ref=ref,
                    component_amount=component_amount,
                ),
            )

        except sqlite3.IntegrityError as e:

            if "UNIQUE constraint failed: component" in str(e):
                mrich.warning(
                    f"Did not add existing component={ref} (type={component_type}) to {route=}"
                )
                return None
            else:
                raise

        component_id = self.cursor.lastrowid

        if commit:
            self.commit()

        return component_id

    def insert_interaction(
        self,
        *,
        feature: Feature | int,
        pose: Pose | int,
        type: str,
        family: str,
        atom_ids: list[int],
        prot_coord: list[float],
        lig_coord: list[float],
        distance: float,
        angle: float | None = None,
        energy: float | None = None,
        warn_duplicate: bool = True,
        commit: bool = True,
        table: str = "interaction",
    ) -> int:
        """Insert an entry into the interaction table

        :param feature: associated :class:`.Feature` object or ID
        :param pose: associated :class:`.Pose` object or ID
        :param type: interaction type
        :param family: ligand feature type
        :param atom_ids: atom indices of ligand feature
        :param prot_coord: ``[x,y,z]`` coordinate of protein feature
        :param lig_coord: ``[x,y,z]`` coordinate of ligand feature
        :param distance: interaction distance ``Angstrom``
        :param angle: optional interaction angle ``degrees``
        :param energy: energy score ``kcal/mol``, defaults to ``None``
        :param warn_duplicate: print a warning if the pose already exists (Default value = True)
        :param commit: commit the changes to the database (Default value = True)
        :param table: the name of the table to insert into (Default value = 'interaction')
        :returns: the interaction ID
        """

        # validation

        if isinstance(feature, Feature):
            feature = feature.id

        if isinstance(pose, Pose):
            pose = pose.id

        from .prolif import FEATURE_FAMILIES, INTERACTION_TYPES

        if family:
            assert family in FEATURE_FAMILIES, f"Unsupported {family=}"
        else:
            family = "Unknown"

        # assert type in INTERACTION_TYPES.values(), f"Unsupported {type=}"

        assert isinstance(atom_ids, list), f"Unsupported {atom_ids=}"
        assert not any(
            [not isinstance(i, int) for i in atom_ids]
        ), f"Unsupported {atom_ids=}"
        atom_ids = json.dumps(atom_ids)

        prot_coord = list(prot_coord) if prot_coord is not None else []
        assert len(prot_coord) == 3 or not prot_coord, f"Unsupported {prot_coord=}"
        assert not any(
            [not isinstance(i, float) for i in prot_coord]
        ), f"Unsupported {prot_coord=}"
        prot_coord = json.dumps(prot_coord)

        lig_coord = list(lig_coord) if lig_coord is not None else []
        assert len(lig_coord) == 3 or not lig_coord, f"Unsupported {lig_coord=}"
        assert not any(
            [not isinstance(i, float) for i in lig_coord]
        ), f"Unsupported {lig_coord=}"
        lig_coord = json.dumps(lig_coord)

        try:
            distance = float(distance)
        except ValueError:
            raise ValueError(f"Unsupported {distance=}")

        try:
            if angle is not None:
                angle = float(angle)
        except ValueError:
            raise ValueError(f"Unsupported {angle=}")

        if energy is not None:
            try:
                energy = float(energy)
            except ValueError:
                raise ValueError(f"Unsupported {energy=}")

        # insertion

        sql = f"""
        INSERT INTO {table}(
            interaction_feature,
            interaction_pose,
            interaction_type,
            interaction_family,
            interaction_atom_ids,
            interaction_prot_coord,
            interaction_lig_coord,
            interaction_distance,
            interaction_angle,
            interaction_energy
        )
        VALUES(?1, ?2, ?3, ?4, ?5, ?6, ?7, ?8, ?9, ?10)
        """

        try:
            self.execute(
                sql,
                (
                    feature,
                    pose,
                    type,
                    family,
                    atom_ids,
                    prot_coord,
                    lig_coord,
                    distance,
                    angle,
                    energy,
                ),
            )

        except sqlite3.IntegrityError as e:
            mrich.error(e)
            if warn_duplicate:
                # mrich.warning(f"Skipping existing interaction: {(feature, pose, family, atom_ids, prot_coord, lig_coord, distance, energy)}")
                mrich.warning(
                    f"Skipping existing interaction: {feature=} {pose=} {family=} {atom_ids=}"
                )
            return None

        except Exception as e:
            mrich.error(e)

        interaction_id = self.cursor.lastrowid

        if commit:
            self.commit()

        return interaction_id

    def insert_subsite(self, target: int, name: str, commit: bool = True) -> int:
        """Insert an entry into the subsite table

        :param target: protein :class:`.Target` ID
        :param name: name of the protein subsite/subsite
        :returns: the subsite ID

        """

        assert isinstance(target, int)
        assert isinstance(name, str)

        sql = """
        INSERT INTO subsite(subsite_target, subsite_name)
        VALUES(?1, ?2)
        """

        try:
            self.execute(sql, (target, name))

        except sqlite3.IntegrityError as e:
            mrich.warning(f"Skipping existing subsite for {target=} with {name=}")
            return None

        except Exception as e:
            mrich.error(e)

        subsite_id = self.cursor.lastrowid
        if commit:
            self.commit()

        return subsite_id

    def insert_subsite_tag(
        self,
        *,
        pose_id: int,
        name: str | None,
        target: int | None = None,
        subsite_id: int | None = None,
        commit: bool = True,
    ) -> int:
        """Insert an entry into the subsite_tag table

        :param pose_id: :class:`.Pose` ID
        :param name: name of the protein subsite/pocket
        :param target: protein :class:`.Target` ID, defaults to querying pose table
        :param target: protein Subsite ID, defaults to querying Subsite table
        :returns: the Subsite ID

        """

        if name is not None:
            assert isinstance(name, str)

        assert isinstance(pose_id, int)

        if not target:
            (target,) = self.select_where(
                table="pose", key="id", value=pose_id, query="pose_target"
            )

        assert isinstance(target, int)

        if not subsite_id:
            subsite_id = self.get_subsite_id(name=name, none="quiet")

        if not subsite_id:
            subsite_id = self.insert_subsite(name=name, target=target)

        assert isinstance(subsite_id, int)

        sql = """
        INSERT INTO subsite_tag(subsite_tag_ref, subsite_tag_pose)
        VALUES(?1, ?2)
        """

        try:
            self.execute(sql, (subsite_id, pose_id))

        except sqlite3.IntegrityError as e:
            mrich.warning(
                f"Skipping existing subsite_tag for {subsite_id=} with {pose_id=}"
            )
            return None

        except Exception as e:
            mrich.error(e)

        subsite_tag_id = self.cursor.lastrowid
        if commit:
            self.commit()

        return subsite_tag_id

    ### SELECTION

    def select(
        self,
        query: str,
        table: str,
        multiple: bool = False,
    ) -> tuple | list[tuple]:
        """Wrapper for the SQL SELECT query, in the following syntax:

        ::

                'SELECT {query} FROM {table}'

        :param query: the columns to return
        :param table: the table from which to select
        :param multiple: fetch all results (Default value = False)
        :returns: the result of the query

        """

        sql = f"SELECT {query} FROM {table}"

        try:
            self.execute(sql)
        except sqlite3.OperationalError as e:
            mrich.var("sql", sql)
            raise

        if multiple:
            result = self.cursor.fetchall()
        else:
            result = self.cursor.fetchone()

        return result

    def select_where(
        self,
        query: str,
        table: str,
        key: str,
        value: str | None = None,
        multiple: bool = False,
        none: str | None = "error",
        sort: str = None,
    ) -> tuple | list[tuple]:
        """Select entries where ``key == value``

        Examples
        ========

        Find compound alias with matching ID:

        ::

                animal.db.select_where(
                        query='compound_alias',
                        table='compound',
                        key='id',
                        value='123',
                )

                # the above evaluates to:
                'SELECT compound_id FROM compound WHERE compound_id = 123'

        Find compound aliases with ID below 10 and order alphabetically:

        ::

                animal.db.select_where(
                        query='compound_alias',
                        table='compound',
                        key='compound_id < 10',
                        multiple=True,
                        sort='compound_alias',
                )

                # the above evaluates to:
                'SELECT compound_id FROM compound WHERE compound_id < 10 ORDER BY compound_alias'

        Parameters
        ==========

        :param query: the columns to return
        :param table: the table from which to select
        :param key: column name to match to value, if no ``value`` is provided the key argument should contain the a SQL string to select entries
        :param value: the value to match (Default value = None)
        :param multiple: fetch all results (Default value = False)
        :param none: define the behaviour for no matches, any value other than ``'error'`` will silently return empty data (Default value = 'error')
        :param sort: optionally sort the output (Default value = None)
        :returns: the result of the query

        """

        if isinstance(value, str):
            if "'" in value:
                value = f'"{value}"'
            else:
                value = f"'{value}'"

        if value is not None:
            where_str = f"{table}_{key}={value}"
        else:
            where_str = key

        if sort:
            sql = f"SELECT {query} FROM {table} WHERE {where_str} ORDER BY {sort}"
        else:
            sql = f"SELECT {query} FROM {table} WHERE {where_str}"

        try:
            self.execute(sql)
        except sqlite3.OperationalError as e:
            mrich.var("sql", sql)
            raise

        if multiple:
            result = self.cursor.fetchall()
        else:
            result = self.cursor.fetchone()

        if not result and none == "error":
            mrich.error(f"No entry in {table} with {where_str}")
            return None
        elif not result and none == "exception":
            raise ValueError(f"No entry in {table} with {where_str}")

        # if not result:
        #     raise ValueError(f"No entry in {table} with {where_str}")

        return result

    def select_id_where(
        self,
        table: str,
        key: str,
        value: str | None = None,
        multiple: bool = False,
        none: str | None = "error",
    ) -> tuple | list[tuple]:
        """Select ID's where ``key==value``. Similar to :meth:`.select_where` except the query argument is always ``{table}_id``.

        :param table: the table from which to select
        :param key: column name to match to value, if no ``value`` is provided this
        :param value: the value to match (Default value = None)
        :param multiple: fetch all results (Default value = False)
        :param none: define the behaviour for no matches, any value other than ``'error'`` will silently return empty data (Default value = 'error')
        :returns: the result of the query

        """
        return self.select_where(
            query=f"{table}_id",
            table=table,
            key=key,
            value=value,
            multiple=multiple,
            none=none,
        )

    def select_all_where(
        self,
        table: str,
        key: str,
        value: str | None = None,
        multiple: bool = False,
        none: str | None = "error",
    ) -> tuple | list[tuple]:
        """Select entries where ``key==value``. Similar to :meth:`.select_where` except the query argument is always ``*``.

        :param table: the table from which to select
        :param key: column name to match to value, if no ``value`` is provided this
        :param value: the value to match (Default value = None)
        :param multiple: fetch all results (Default value = False)
        :param none: define the behaviour for no matches, any value other than ``'error'`` will silently return empty data (Default value = 'error')
        :returns: the result of the query

        """
        return self.select_where(
            query="*", table=table, key=key, value=value, multiple=multiple, none=none
        )

    ### DELETION

    def delete_where(
        self,
        table: str,
        key: str,
        value: str | None = None,
        commit: bool = True,
    ) -> None:
        """Delete entries where ``key==value``

        :param table: the table from which to delete
        :param key: column name to match to value, if no ``value`` is provided this
        :param value: the value to match (Default value = None)
        :param commit: commit the changes (Default value = True)

        """

        if value is not None:

            if isinstance(value, str):
                value = f"'{value}'"

            sql = f"DELETE FROM {table} WHERE {table}_{key}={value}"

        else:

            sql = f"DELETE FROM {table} WHERE {key}"

        try:
            result = self.execute(sql)

        except sqlite3.OperationalError as e:
            mrich.var("sql", sql)
            raise

        if commit:
            self.commit()

    def delete_tag(
        self,
        tag: str,
    ) -> None:
        """Delete all tag entries with the matching name

        :param tag: tag name to match

        """
        self.delete_where(table="tag", key="name", value=tag)

    def delete_interactions(self) -> None:
        """Delete all calculated interactions and set pose_fingerprint appropriately"""

        self.delete_where(table="interaction", key="interaction_id > 0")
        self.update_all(table="pose", key="pose_fingerprint", value=0)

    def delete_features(self) -> None:
        """Delete all protein features"""

        self.delete_where(table="feature", key="feature_id > 0")

    def delete_reactions(self) -> None:
        """Delete all reaction data"""

        tables = ["reaction", "reactant", "route", "component"]

        for table in tables:
            self.execute(f"DELETE FROM {table};")
        self.commit()

    def delete_subsites(self) -> None:
        """Delete all protein subsites"""

        self.delete_where(table="subsite", key="subsite_id > 0")
        self.delete_where(table="subsite_tag", key="subsite_tag_id > 0")

    ### UPDATE

    def update(
        self,
        *,
        table: str,
        id: int,
        key: str,
        value,
        commit: bool = True,
    ) -> int:
        """Update a field in a database entry with given ID

        :param table: the table which to update
        :param id: the ID of the entry to update
        :param key: column name to update
        :param value: the value to insert
        :param commit: commit the changes to the database (Default value = True)
        :returns: the ID of the modified entry

        """

        sql = f"""
        UPDATE {table}
        SET {key} = ?
        WHERE {table}_id = {id};
        """

        try:
            self.execute(sql, (value,))
        except sqlite3.OperationalError as e:
            mrich.var("sql", sql)
            raise

        id = self.cursor.lastrowid
        if commit:
            self.commit()
        return id

    def update_all(
        self,
        *,
        table: str,
        key: str,
        value,
        commit: bool = True,
    ) -> None:
        """Update all fields in a table column at once

        :param table: the table which to update
        :param key: column name to update
        :param value: the value to insert
        :param commit: commit the changes to the database (Default value = True)

        """

        sql = f"""
        UPDATE {table}
        SET {key} = ?
        """

        try:
            self.execute(sql, (value,))
        except sqlite3.OperationalError as e:
            mrich.var("sql", sql)
            raise

        if commit:
            self.commit()

    ### COPYING / MIGRATION

    def copy_temp_interactions(self) -> int:
        """Copy the records from the 'temp_interaction' table to the 'interaction' table

        :returns: ID of the last inserted :class:`.Interaction`
        """

        cursor = self.execute(
            """
            INSERT OR IGNORE INTO interaction(
                interaction_feature, 
                interaction_pose, 
                interaction_type, 
                interaction_family, 
                interaction_atom_ids, 
                interaction_prot_coord, 
                interaction_lig_coord, 
                interaction_distance, 
                interaction_angle, 
                interaction_energy
            )
            SELECT interaction_feature, 
                interaction_pose, 
                interaction_type, 
                interaction_family, 
                interaction_atom_ids, 
                interaction_prot_coord, 
                interaction_lig_coord, 
                interaction_distance, 
                interaction_angle, 
                interaction_energy 
            FROM temp_interaction
        """
        )

        return cursor.lastrowid

    def copy_interactions_to_temp(self, pose_id: int) -> int:
        """Copy the records from the 'temp_interaction' table to the 'interaction' table

        :returns: ID of the last inserted :class:`.Interaction`
        """

        cursor = self.execute(
            f"""
            INSERT OR IGNORE INTO temp_interaction(
                interaction_feature, 
                interaction_pose, 
                interaction_type, 
                interaction_family, 
                interaction_atom_ids, 
                interaction_prot_coord, 
                interaction_lig_coord, 
                interaction_distance, 
                interaction_angle, 
                interaction_energy
            )
            SELECT interaction_feature, 
                interaction_pose, 
                interaction_type, 
                interaction_family, 
                interaction_atom_ids, 
                interaction_prot_coord, 
                interaction_lig_coord, 
                interaction_distance, 
                interaction_angle, 
                interaction_energy 
            FROM interaction
            WHERE interaction_pose = {pose_id}
        """
        )

        return cursor.lastrowid

    def migrate_legacy_bases(self) -> int:
        """Migrate legacy compound_base records from the 'compound' table to the 'scaffold' table

        :returns: ID of the last inserted scaffold record
        """

        mrich.debug("HIPPO.Database.migrate_legacy_bases()")

        cursor = self.execute(
            """
            INSERT INTO scaffold(scaffold_base, scaffold_superstructure)
            SELECT compound_base, compound_id FROM compound
            WHERE compound_base IS NOT NULL
            """
        )

        self.commit()

        return cursor.lastrowid

    def update_legacy_routes(self) -> None:
        """Update legacy component entries"""

        # add column

        sql = """
        ALTER TABLE component
        ADD component_amount REAL;
        """

        self.execute(sql)

        # set values

        sql = """
        UPDATE component
        SET component_amount = :component_amount
        WHERE component_type = :component_type;
        """

        self.execute(sql, dict(component_amount=None, component_type=1))
        self.execute(sql, dict(component_amount=1.0, component_type=2))
        self.execute(sql, dict(component_amount=1.0, component_type=3))

    def update_legacy_reaction_metadata(self) -> None:
        """Add reaction_metadata column"""

        sql = """
        ALTER TABLE reaction
        ADD reaction_metadata TEXT;
        """

        self.execute(sql)

    def update_compound_pattern_bfp_table(self):
        """Update the compound pattern BFP table"""
        self.execute(
            """
            INSERT INTO compound_pattern_bfp
            SELECT c.compound_id, c.compound_pattern_bfp FROM compound AS c
            LEFT JOIN compound_pattern_bfp as fp
            ON c.compound_id = fp.compound_id
            WHERE fp.compound_id IS NULL
        """
        )

    def prune_duplicate_routes(self) -> None:
        """Remove duplicate routes from the database"""

        from collections import Counter

        sql = """
        SELECT route_id, route_product, component_ref, component_type FROM route
        INNER JOIN component ON route_id = component_route
        """

        records = self.execute(sql).fetchall()

        routes = {}

        for route_id, route_product, component_ref, component_type in records:
            if route_id not in routes:
                routes[route_id] = (route_product, set())

            key = (component_ref, component_type)
            routes[route_id][1].add(key)

        routes = {key: (value[0], tuple(value[1])) for key, value in routes.items()}

        flat_routes = [value for key, value in routes.items()]

        mrich.var("#routes", len(flat_routes))

        counter = Counter(flat_routes)
        duplicates = {item: count for item, count in counter.items() if count > 1}

        mrich.var("products with duplicate routes", len(duplicates))

        if not duplicates:
            mrich.success("No duplicate routes found")
            return None

        delete = set()
        for dupe in duplicates:
            matched_ids = [k for k, v in routes.items() if v == dupe]
            mrich.print(
                "compound", dupe[0], "has", len(matched_ids), "duplicate routes"
            )
            for route_id in matched_ids[1:]:
                delete.add(route_id)

        str_ids = str(tuple(delete)).replace(",)", ")")
        self.delete_where(table="component", key=f"component_route IN {str_ids}")
        self.delete_where(table="route", key=f"route_id IN {str_ids}")

        mrich.success("Deleted", len(delete), "duplicate routes")

        return delete

    def calculate_all_scaffolds(self):
        """Determine and insert records for all substructure/superstructure relationships in the Compound table"""

        n_before = self.count("scaffold")

        mrich.var("#compounds", self.count("compound"))
        mrich.var("#scaffold defs", n_before)

        sql = """
        SELECT compound_id, compound_mol, compound_pattern_bfp 
        FROM compound
        """

        with mrich.loading("Fetching compounds..."):
            records = self.execute(sql).fetchall()

        sql = """
            INSERT OR IGNORE INTO scaffold
            SELECT ?1, c.compound_id 
            FROM compound AS c, compound_pattern_bfp AS fp
            WHERE c.compound_id = fp.compound_id
            AND c.compound_id <> ?1
            AND mol_is_substruct(c.compound_mol, ?2)
            AND fp.compound_id MATCH rdtree_subset(?3)
        """

        with mrich.loading("Calculating scaffolds..."):
            t1 = time.time()
            self.executemany(sql, records)
            mrich.print("Took", f"{time.time() - t1:.1f}", "seconds")

        with mrich.loading("Committing..."):
            self.commit()

        diff = self.count("scaffold") - n_before

        if diff:
            mrich.success(
                "Found", diff, "new substructure-superstructure relationships"
            )
        else:
            mrich.warning(
                "Found", diff, "new substructure-superstructure relationships"
            )

    def reinitialise_molecules(self):
        """In the case where the Mol binaries in a database are throwing unpickling errors, run this to reinitialise them all from their smiles."""

        mrich.var("#compounds", self.count("compound"))

        sql = """
        UPDATE your_table_name
        SET compound_mol = mol_from_smiles(compound_smiles);
        """

        with mrich.loading("Reinitialising compounds..."):
            self.execute(sql)

        mrich.success("compound_mol records updated")

    def register_compounds(
        self,
        *,
        smiles: list[str],
        radical: str = "warning",
        sanitisation_verbosity: bool = True,
        sanitise: bool = True,
        debug: bool = False,
    ) -> list[tuple[str, str]]:
        """Bulk register compounds"""

        values = []

        if len(smiles) > 1000:
            generator = mrich.track(smiles, prefix="Sanitising...")
        else:
            generator = smiles

        for s in generator:

            if sanitise:
                try:
                    new_smiles = sanitise_smiles(
                        s,
                        sanitisation_failed="error",
                        radical=radical,
                        verbosity=sanitisation_verbosity,
                    )
                except SanitisationError as e:
                    mrich.error(f"Could not sanitise {s=}")
                    mrich.error(str(e))
                    continue
                except AssertionError:
                    mrich.error(f"Could not sanitise {s=}")
                    continue
            else:
                new_smiles = s

            inchikey = inchikey_from_smiles(new_smiles)
            values.append((inchikey, new_smiles))

        sql = """
        INSERT OR IGNORE INTO compound(compound_inchikey, compound_smiles, compound_mol, compound_pattern_bfp, compound_morgan_bfp)
        VALUES(?1, ?2, mol_from_smiles(?2), mol_pattern_bfp(mol_from_smiles(?2), 2048), mol_morgan_bfp(mol_from_smiles(?2), 2, 2048))
        """

        if debug:
            mrich.debug("Inserting...")

        self.executemany(sql, values)
        self.commit()

        self.update_compound_pattern_bfp_table()

        return values

    def register_poses(self, dicts: list[dict]) -> set[int]:
        """Insert or ignore a bunch of poses, also returns a set of Pose IDs"""

        from json import dumps

        ### POSES

        sql = """
        INSERT OR IGNORE INTO pose(
            pose_inchikey, 
            pose_smiles, 
            pose_alias, 
            pose_reference, 
            pose_path, 
            pose_compound, 
            pose_target, 
            pose_mol, 
            pose_energy_score, 
            pose_distance_score, 
            pose_metadata
        )
        VALUES(?1, ?2, ?3, ?4, ?5, ?6, ?7, ?8, ?9, ?10, ?11)
        """

        values = []
        for i, d in enumerate(dicts):

            alias = d["alias"]

            if not alias:
                alias = None

            try:
                values.append(
                    (
                        str(d["inchikey"]),
                        str(d["smiles"]),
                        alias,
                        int(d["reference_id"]),
                        str(d["path"]),
                        int(d["compound_id"]),
                        int(d["target_id"]),
                        d["mol"].ToBinary(),
                        float(d["energy_score"]),
                        float(d["distance_score"]),
                        dumps(d["metadata"]),
                    )
                )
            except KeyError as e:
                mrich.error("Skipping", i, str(e))

        self.executemany(sql, values)
        self.commit()

        ### INSPIRATIONS

        lookup = self.get_pose_path_id_dict()

        values = []
        pose_ids = set()

        for i, d in enumerate(dicts):
            if "inspiration_ids" not in d:
                continue
            derivative_id = lookup.get(str(d["path"]))
            if not derivative_id:
                continue
            pose_ids.add(derivative_id)
            for inspiration_id in d["inspiration_ids"]:
                values.append((inspiration_id, derivative_id))

        sql = """
        INSERT OR IGNORE INTO inspiration(inspiration_original, inspiration_derivative)
        VALUES(?1, ?2)
        """

        self.executemany(sql, values)
        self.commit()

        return pose_ids

    def calculate_all_murcko_scaffolds(self, generic: bool = True):
        """Determine Murcko and optionally generic Murcko scaffolds for all Compounds in the Database and add relevant records.

        :param generic: Calculate generic (single bonds and all carbon) scaffolds as well
        """

        n_before = self.count("scaffold")

        mrich.var("#compounds", self.count("compound"))
        mrich.var("#scaffold defs", n_before)

        from rdkit.Chem import MolFromSmiles, MolToSmiles
        from rdkit.Chem.Scaffolds.MurckoScaffold import (
            MurckoScaffoldSmiles,
            MakeScaffoldGeneric,
        )

        compound_records = self.select(
            query="compound_id, compound_smiles", table="compound", multiple=True
        )

        ### CALCULATE SCAFFOLDS

        murcko_data = {}
        generic_data = {}
        generic_to_murcko = {}
        for c_id, smiles in mrich.track(compound_records):

            # murcko

            try:
                murcko_smiles = sanitise_smiles(MurckoScaffoldSmiles(smiles))
            except KeyboardInterrupt:
                raise
            except:
                mrich.error("can't make murcko:", smiles)
                continue

            if murcko_smiles not in murcko_data:
                murcko_data[murcko_smiles] = set()

            murcko_data[murcko_smiles].add(c_id)

            # generic

            if generic:

                try:
                    generic_smiles = sanitise_smiles(
                        MolToSmiles(MakeScaffoldGeneric(MolFromSmiles(murcko_smiles)))
                    )
                except KeyboardInterrupt:
                    raise
                except:
                    mrich.error("can't make generic:", murcko_smiles)
                    continue

                if generic_smiles not in generic_data:
                    generic_data[generic_smiles] = set()

                if generic_smiles not in generic_to_murcko:
                    generic_to_murcko[generic_smiles] = set()

                generic_data[generic_smiles].add(c_id)
                generic_to_murcko[generic_smiles].add(murcko_smiles)

        mrich.var("#murcko scaffolds", len(murcko_data))
        mrich.var("#generic murcko scaffolds", len(generic_data))

        ### REGISTER MURCKOS

        murcko_values = self.register_compounds(
            smiles=murcko_data.keys(), sanitisation_verbosity=False, sanitise=False
        )
        murcko_s2i = {s: i for i, s in murcko_values}

        ### REGISTER GENERICS

        if generic:
            generic_values = self.register_compounds(
                smiles=generic_data.keys(), sanitisation_verbosity=False, sanitise=False
            )
            generic_s2i = {s: i for i, s in generic_values}

        ### TAG MURCKOS

        murcko_ids = self.select_id_where(
            table="compound",
            key=f"compound_inchikey IN {tuple(murcko_s2i.values())}",
            multiple=True,
        )

        self.executemany(
            """INSERT OR IGNORE INTO tag (tag_name, tag_compound) VALUES (?,?)""",
            [("MurckoScaffold", i) for i, in murcko_ids],
        )

        ### TAG GENERICS

        if generic:

            generic_ids = self.select_id_where(
                table="compound",
                key=f"compound_inchikey IN {tuple(generic_s2i.values())}",
                multiple=True,
            )

            self.executemany(
                """INSERT OR IGNORE INTO tag (tag_name, tag_compound) VALUES (?,?)""",
                [("GenericMurckoScaffold", i) for i, in generic_ids],
            )

        ### ADD MURCKO SCAFFOLD RELATIONS

        pairs = []

        murcko_inchikey_lookup = self.get_compound_inchikey_id_dict(murcko_s2i.values())
        for murcko_smiles, c_ids in murcko_data.items():
            murcko_id = murcko_inchikey_lookup[murcko_s2i[murcko_smiles]]
            for c_id in c_ids:
                pairs.append((murcko_id, c_id))

        pairs = [(a, b) for a, b in pairs if a != b]

        mrich.var("#murcko scaffold relations", len(pairs))

        self.executemany(
            """INSERT OR IGNORE INTO scaffold (scaffold_base, scaffold_superstructure) VALUES (?,?)""",
            pairs,
        )

        ### ADD GENERIC SCAFFOLD RELATIONS

        if generic:

            pairs = []

            generic_inchikey_lookup = self.get_compound_inchikey_id_dict(
                generic_s2i.values()
            )
            for generic_smiles, c_ids in generic_data.items():
                generic_id = generic_inchikey_lookup[generic_s2i[generic_smiles]]
                for c_id in c_ids:
                    pairs.append((generic_id, c_id))

            for generic_smiles, murcko_smiles_list in generic_to_murcko.items():
                generic_id = generic_inchikey_lookup[generic_s2i[generic_smiles]]
                for murcko_smiles in murcko_smiles_list:
                    murcko_id = murcko_inchikey_lookup[murcko_s2i[murcko_smiles]]
                    pairs.append((generic_id, murcko_id))

            pairs = [(a, b) for a, b in pairs if a != b]

            mrich.var("#generic murcko scaffold relations", len(pairs))

            self.executemany(
                """INSERT OR IGNORE INTO scaffold (scaffold_base, scaffold_superstructure) VALUES (?,?)""",
                pairs,
            )

        if generic:
            return murcko_data, generic_data
        else:
            return murcko_data

    def set_derivative_subsites(self, commit: bool = True) -> None:
        """Propagate all subsite assignments from inspirations to their derivatives"""

        sql = """
        INSERT OR IGNORE INTO subsite_tag(subsite_tag_ref, subsite_tag_pose)
        SELECT subsite_tag_ref, inspiration_derivative FROM subsite_tag
        INNER JOIN inspiration
        ON subsite_tag_pose = inspiration_original
        """

        self.execute(sql)

        if commit:
            self.commit()

    ### GETTERS

    def get_compound(
        self,
        *,
        id: int | None = None,
        inchikey: str | None = None,
        alias: str | None = None,
        smiles: str | None = None,
        none: str = "error",
        **kwargs,
    ) -> Compound:
        """Get a :class:`.Compound` using one of the following fields: ['id', 'inchikey', 'alias', 'smiles']

        :param id: the ID to search for (Default value = None)
        :param inchikey: the InChi-Key to search for (Default value = None)
        :param alias: the alias to search for (Default value = None)
        :param smiles: the smiles to search for (Default value = None)
        :returns: the :class:`.Compound` object

        """

        if id is None:
            id = self.get_compound_id(
                inchikey=inchikey, smiles=smiles, alias=alias, none=none, **kwargs
            )

        if not id:
            if none == "error":
                mrich.error(f"Invalid {id=}")
            return None

        query = "compound_id, compound_inchikey, compound_alias, compound_smiles"
        entry = self.select_where(
            query=query, table="compound", key="id", value=id, none=none, **kwargs
        )
        compound = Compound(self._animal, self, *entry, metadata=None, mol=None)
        return compound

    def get_compound_id(
        self,
        *,
        inchikey: str | None = None,
        alias: str | None = None,
        smiles: str | None = None,
        **kwargs,
    ) -> int:
        """Get a compound's ID using one of the following fields: ['inchikey', 'alias', 'smiles']

        :param inchikey: the InChi-Key to search for (Default value = None)
        :param alias: the alias to search for (Default value = None)
        :param smiles: the smiles to search for (Default value = None)
        :returns: the :class:`.Compound` ID

        """

        if inchikey:
            entry = self.select_id_where(
                table="compound", key="inchikey", value=inchikey, **kwargs
            )

        elif alias:
            entry = self.select_id_where(
                table="compound", key="alias", value=alias, **kwargs
            )

        elif smiles:
            entry = self.select_id_where(
                table="compound", key="smiles", value=smiles, **kwargs
            )

        else:
            raise NotImplementedError

        if entry:
            return entry[0]

        return None

    def get_compound_computed_property(
        self,
        prop: str,
        compound_id: int,
    ) -> int | str:
        """Use chemicalite to calculate a property from the stored binary molecule

        :param prop: the property to calculate [num_heavy_atoms, formula, num_rings]
        :param compound_id: the compound ID to query
        :returns: the value of the computed property

        """

        function = CHEMICALITE_COMPOUND_PROPERTY_MAP[prop]
        (val,) = self.select_where(
            query=f"{function}(compound_mol)",
            table="compound",
            key="id",
            value=compound_id,
            multiple=False,
        )
        return val

    def get_pose(
        self,
        *,
        id: int | None = None,
        inchikey: str = None,
        alias: str = None,
    ) -> Pose:
        """Get a pose using one of the following fields: ['id', 'inchikey', 'alias']

        :param id: the ID to search for (Default value = None)
        :param inchikey: the InChi-Key to search for (Default value = None)
        :param alias: the alias to search for (Default value = None)
        :returns: the :class:`.Pose` object

        """

        if id is None:
            id = self.get_pose_id(inchikey=inchikey, alias=alias)

        if isinstance(id, list):
            from .pset import PoseSet

            return PoseSet(self, id)

        if not id:
            mrich.error(f"Invalid {id=}")
            return None

        query = "pose_id, pose_inchikey, pose_alias, pose_smiles, pose_reference, pose_path, pose_compound, pose_target, pose_mol, pose_fingerprint, pose_energy_score, pose_distance_score"
        entry = self.select_where(query=query, table="pose", key="id", value=id)
        pose = Pose(self, *entry)
        return pose

    def get_pose_id(
        self,
        *,
        inchikey: str | None = None,
        alias: str | None = None,
    ) -> int:
        """Get a pose's ID using one of the following fields: ['inchikey', 'alias', 'smiles']

        :param table: the table from which to get the entry (Default value = 'pose')
        :param inchikey: the InChi-Key to search for (Default value = None)
        :param alias: the alias to search for (Default value = None)
        :returns: the :class:`.Pose` ID

        """

        if inchikey:
            # inchikey might not be unique
            entries = self.select_id_where(
                table="pose", key="inchikey", value=inchikey, multiple=True
            )
            if len(entries) != 1:
                mrich.warning(f"Multiple poses with {inchikey=}")
                return [i for i, in entries]
            else:
                entry = entries[0]

        elif alias:
            entry = self.select_id_where(table="pose", key="alias", value=alias)

        else:
            raise NotImplementedError

        if entry:
            return entry[0]

        return None

    def get_reaction(
        self,
        *,
        id: int | None = None,
        none: str | None = None,
    ) -> Reaction:
        """Get a reaction using its ID

        :param id: the ID to search for (Default value = None)
        :param none: define the behaviour for no matches, any value other than ``'error'`` will silently return empty data (Default value = 'error')
        :returns: the :class:`.Reaction` object

        """

        if not id:
            mrich.error(f"Invalid {id=}")
            return None

        query = "reaction_id, reaction_type, reaction_product, reaction_product_yield"
        entry = self.select_where(
            query=query, table="reaction", key="id", value=id, none=none
        )

        if not entry:
            return None

        reaction = Reaction(self, *entry)
        return reaction

    def get_quote(
        self,
        *,
        id: int | None = None,
        none: str | None = None,
    ) -> Quote:
        """Get a quote using its ID

        :param id: the ID to search for (Default value = None)
        :param none: define the behaviour for no matches, any value other than ``'error'`` will silently return empty data (Default value = 'error')
        :returns: the :class:`.Quote` object

        """

        query = "quote_compound, quote_supplier, quote_catalogue, quote_entry, quote_amount, quote_price, quote_currency, quote_lead_time, quote_purity, quote_date, quote_smiles, quote_id "
        entry = self.select_where(
            query=query, table="quote", key="id", value=id, none=none
        )

        return Quote(
            db=self,
            id=entry[11],
            compound=entry[0],
            supplier=entry[1],
            catalogue=entry[2],
            entry=entry[3],
            amount=entry[4],
            price=entry[5],
            currency=entry[6],
            lead_time=entry[7],
            purity=entry[8],
            date=entry[9],
            smiles=entry[10],
        )

    def get_metadata(
        self,
        *,
        table: str,
        id: int,
    ) -> dict:
        """Get metadata dictionary from a specific table and ID

        :param table: the table from which to get the entry
        :param id: the ID to search for (Default value = None)
        :returns: a dictionary of metadata

        """

        (payload,) = self.select_where(
            query=f"{table}_metadata", table=table, key=f"id", value=id
        )

        if payload:
            payload = json.loads(payload)

        else:
            payload = dict()

        metadata = MetaData(payload)

        metadata._db = self
        metadata._id = id
        metadata._table = table

        return metadata

    def get_target(
        self,
        *,
        id: int,
    ) -> Target:
        """Get target with specific ID

        :param id: the ID of the target to retrieve
        :returns: :class:`.Target` object

        """

        return Target(db=self, id=id, name=self.get_target_name(id=id))

    def get_target_name(
        self,
        *,
        id: int,
    ) -> str:
        """Get the name of a target with given ID

        :param id: the ID of the target to retrieve
        :returns: target name

        """

        table = "target"
        (payload,) = self.select_where(
            query=f"{table}_name", table=table, key="id", value=id
        )
        return payload

    def get_target_id(
        self,
        *,
        name: str,
    ) -> int | None:
        """Get target ID with a given name

        :param name: the protein target name
        :returns: the :class:`.Target` ID

        """

        table = "target"
        entry = self.select_id_where(table=table, key="name", value=name)

        if entry:
            return entry[0]

        return None

    def get_feature(
        self,
        *,
        id: int,
    ) -> Feature:
        """Get a protein interaction :class:`.Feature` with a given ID

        :param id: the protein interaction :class:`.Feature` ID to be retrieved
        :returns: :class:`.Feature` object

        """

        entry = self.select_all_where(table="feature", key="id", value=id)
        return Feature(*entry)

    def get_route(
        self,
        *,
        id: int,
        debug: bool = False,
    ) -> Route:
        """Fetch a :class:`.Route` object stored in the :class:`.Database`.

        :param id: the ID of the :class:`.Route` to be retrieved
        :param debug: increase verbosity for debugging, defaults to False
        :returns: :class:`.Route` object

        """

        from .cset import CompoundSet, IngredientSet
        from .rset import ReactionSet

        (product_id,) = self.select_where(
            table="route", query="route_product", key="id", value=id
        )

        if debug:
            mrich.var("product_id", product_id)

        triples = self.select_where(
            table="component",
            query="component_ref, component_type, component_amount",
            key=f"component_route IS {id} ORDER BY component_id",
            multiple=True,
        )

        reaction_ids = []
        reactant_ids = []
        reactant_amounts = []
        intermediate_ids = []
        intermediate_amounts = []

        for ref, c_type, amount in triples:
            match c_type:
                case 1:
                    reaction_ids.append(ref)
                case 2:
                    reactant_ids.append(ref)
                    reactant_amounts.append(amount)
                case 3:
                    intermediate_ids.append(ref)
                    intermediate_amounts.append(amount)
                case _:
                    raise ValueError(f"Unknown component type {c_type}")

        if debug:
            mrich.var("pairs", pairs)

        products = CompoundSet(self, [product_id])
        reactants = CompoundSet(self, reactant_ids)
        intermediates = CompoundSet(self, intermediate_ids)

        products = IngredientSet.from_compounds(compounds=products, amount=1)
        reactants = IngredientSet.from_compounds(
            compounds=reactants, amount=reactant_amounts
        )
        intermediates = IngredientSet.from_compounds(
            compounds=intermediates, amount=intermediate_amounts
        )

        reactions = ReactionSet(self, reaction_ids)

        recipe = Route(
            self,
            route_id=id,
            product=products,
            reactants=reactants,
            intermediates=intermediates,
            reactions=reactions,
        )

        if debug:
            mrich.var("recipe", recipe)

        return recipe

    def get_route_products(self) -> "CompoundSet | None":
        """Get a :class:`.CompoundSet` of all route products"""
        from .cset import CompoundSet

        records = self.execute("SELECT DISTINCT route_product FROM route").fetchall()
        if not records:
            return None
        return CompoundSet(self, [i for i, in records])

    def get_route_id_product_dict(self) -> dict[int, int]:
        """Get a dictionary mapping route ID's to their product :class:`.Compound`"""
        records = self.execute("SELECT route_id, route_product FROM route").fetchall()
        return {route_id: route_product for route_id, route_product in records}

    def get_product_id_routes_dict(self) -> dict[int, set[int]]:
        """Get a dictionary mapping product :class:`.Compound` to their route IDs"""
        records = self.execute("SELECT route_id, route_product FROM route").fetchall()

        lookup = {}

        for route_id, route_product in records:

            if route_product not in lookup:
                lookup[route_product] = set()

            lookup[route_product].add(route_id)

        return lookup

    def get_compound_id_pose_ids_dict(self, cset: "CompoundSet") -> dict[int, set]:
        """Get a dictionary mapping :class:`.Compound` ID's to their associated :class:`.Pose` ID's"""
        records = self.execute(
            f"SELECT pose_compound, pose_id FROM pose WHERE pose_compound IN {cset.str_ids}"
        ).fetchall()

        d = {}

        for comp_id, pose_id in records:
            d[comp_id] = d.get(comp_id, set())
            d[comp_id].add(pose_id)
        return d

    def get_compound_id_suppliers_dict(
        self, cset: "CompoundSet"
    ) -> dict[int, set[str]]:
        """Get a dictionary mapping :class:`.Compound` ID's to suppliers which stock it"""
        records = self.execute(
            f"SELECT quote_compound, quote_supplier FROM quote WHERE quote_compound IN {cset.str_ids}"
        ).fetchall()

        d = {}

        for comp_id, quote_supplier in records:
            d[comp_id] = d.get(comp_id, set())
            d[comp_id].add(quote_supplier)

        return d

    def get_compound_id_smiles_dict(self, cset: "CompoundSet") -> dict[int, set[str]]:
        """Get a dictionary mapping :class:`.Compound` ID's to suppliers which stock it"""
        records = self.execute(
            f"SELECT compound_id, compound_smiles FROM compound WHERE compound_id IN {cset.str_ids}"
        ).fetchall()

        d = {}

        for comp_id, comp_smiles in records:
            d[comp_id] = comp_smiles

        return d

    def get_compound_inchikey_id_dict(self, inchikeys: list[str]):
        """Get a dictionary mapping :class:`.Compound` inchikeys to their ID's"""

        inchikey_str = str(tuple(inchikeys)).replace(",)", ")")

        records = self.select_where(
            table="compound",
            multiple=True,
            query="compound_inchikey, compound_id",
            key=f"compound_inchikey IN {inchikey_str}",
        )

        return {
            compound_inchikey: compound_id for compound_inchikey, compound_id in records
        }

    def get_compound_cluster_dict(
        self,
        cset: "CompoundSet | None" = None,
        *,
        fractions: bool = False,
        max_bases: int | None = None,
        fraction_reference: "CompoundSet | None" = None,
    ) -> dict[tuple, set]:
        """Create a dictionary grouping compounds by their scaffold/base cluster.

        :param cset: :class:`.CompoundSet` subset to query, defaults to all compounds
        :param fractions: Calculate fractional populations for each cluster
        :param max_bases: Define the maximum number of compounds to use as cluster keys
        :param fraction_reference: Use cset to build the cluster map and use fraction_reference to determine the fractional populations
        :returns: A dictionary mapping a tuple of scaffold :class:`.Compound` IDs to a set of superstructure :class:`.Compound` ID's.
        """

        if fractions:
            assert cset is not None

        if (cset is not None and not fractions) or (fraction_reference is not None):
            sql = f"""
            SELECT scaffold_superstructure, scaffold_base FROM scaffold
            WHERE scaffold_superstructure IN {cset.str_ids}
            """
        else:
            sql = f"""
            SELECT scaffold_superstructure, scaffold_base FROM scaffold
            """

        records = self.execute(sql).fetchall()

        lookup = {}
        for superstructure, base in records:
            group = lookup.setdefault(superstructure, set())
            group.add(base)

        clustered = {}
        for superstructure, bases in lookup.items():
            if max_bases and len(bases) > max_bases:
                continue
            cluster = tuple(bases)
            group = clustered.setdefault(cluster, set())
            group.add(superstructure)

        if fractions:
            if fraction_reference is not None:
                cset = fraction_reference

            fractions = {}
            for cluster, members in clustered.items():
                size = len(members)
                present = sum(1 for c in members if c in cset)
                fractions[cluster] = present / size
            return fractions

        return clustered

    def get_compound_scaffold_dict(self) -> dict[int, set[int]]:
        """Get a dictionary mapping scaffold_base compound ID's to a set of their superstructure IDs"""

        records = self.select(
            table="scaffold",
            query="scaffold_base, scaffold_superstructure",
            multiple=True,
        )

        data = {}
        for base, elab in records:
            if base not in data:
                data[base] = set()
            data[base].add(elab)

        return data

    def get_compound_tag_dict(
        self,
        cset: "CompoundSet | None" = None,
    ) -> dict[int, set[str]]:
        """Get a dictionary mapping compound ID's to their tags"""

        if cset:
            raise NotImplementedError

        records = self.select(
            query="tag_name, tag_compound", table="tag", multiple=True
        )

        data = {}
        for tag_name, compound_id in records:
            if compound_id not in data:
                data[compound_id] = set()
            data[compound_id].add(tag_name)

        # null IDS
        comp_ids = self.select(table="compound", query="compound_id", multiple=True)
        comp_ids = set(q for q, in comp_ids)

        null_ids = comp_ids - set(data.keys())

        for c_id in null_ids:
            data[c_id] = set()

        return data

    def get_pose_tag_dict(
        self,
        pset: "PoseSet | None" = None,
    ) -> dict[int, set[str]]:
        """Get a dictionary mapping pose ID's to their tags"""

        if pset:
            raise NotImplementedError

        records = self.select(query="tag_name, tag_pose", table="tag", multiple=True)

        data = {}
        for tag_name, pose_id in records:
            if pose_id not in data:
                data[pose_id] = set()
            data[pose_id].add(tag_name)

        # null IDS
        pose_ids = self.select(table="pose", query="pose_id", multiple=True)
        pose_ids = set(q for q, in pose_ids)

        null_ids = pose_ids - set(data.keys())

        for c_id in null_ids:
            data[c_id] = set()

        return data

    def get_pose_subsite_names_dict(self) -> dict[int, set[str]]:
        """Get a dictionary mapping pose ID's to their subsite names"""

        lookup = {
            i: n
            for i, n in self.select(
                query="subsite_id, subsite_name", table="subsite", multiple=True
            )
        }

        records = self.select(
            query="subsite_tag_ref, subsite_tag_pose",
            table="subsite_tag",
            multiple=True,
        )

        data = {}
        for subsite_id, pose_id in records:
            data.setdefault(pose_id, set())
            data[pose_id].add(lookup[subsite_id])

        # null IDS
        pose_ids = self.select(table="pose", query="pose_id", multiple=True)
        pose_ids = set(q for q, in pose_ids)

        null_ids = pose_ids - set(data.keys())

        for c_id in null_ids:
            data[c_id] = set()

        return data

    def get_pose_id_interaction_ids_dict(self, pset: "PoseSet") -> dict[int, set]:
        """Get a dictionary mapping :class:`.Pose` ID's to their associated :class:`.Interaction` ID's"""
        records = self.execute(
            f"SELECT interaction_pose, interaction_id FROM interaction WHERE interaction_pose IN {pset.str_ids}"
        ).fetchall()

        d = {}

        for pose_id, interaction_id in records:
            d[pose_id] = d.get(pose_id, set())
            d[pose_id].add(interaction_id)
        return d

    def get_pose_alias_id_dict(self, pset: "PoseSet | None" = None) -> dict[str, int]:
        """Get a dictionary mapping :class:`.Pose` aliases to ID's"""

        if pset:
            records = self.execute(
                f"""
            SELECT pose_id, pose_alias FROM pose 
            WHERE pose_alias IS NOT NULL
            AND pose_id IN {pset.str_ids}"""
            ).fetchall()

        else:
            records = self.execute(
                """SELECT pose_id, pose_alias FROM pose 
            WHERE pose_alias IS NOT NULL"""
            ).fetchall()

        d = {}
        for pose_id, pose_alias in records:
            d[pose_alias] = pose_id

        return d

    def get_pose_alias_path_dict(self, pset: "PoseSet | None" = None) -> dict[str, str]:
        """Get a dictionary mapping :class:`.Pose` aliases to paths"""

        if pset:
            records = self.execute(
                f"""
            SELECT pose_alias, pose_path FROM pose 
            AND pose_id IN {pset.str_ids}"""
            ).fetchall()

        else:
            records = self.execute(
                """SELECT pose_alias, pose_path FROM pose """
            ).fetchall()

        d = {}
        for pose_alias, pose_path in records:
            d[pose_alias] = pose_path

        return d

    def get_pose_id_alias_dict(self, pset: "PoseSet | None" = None) -> dict[str, int]:
        """Get a dictionary mapping :class:`.Pose` aliases to ID's"""

        if pset:
            records = self.execute(
                f"""
            SELECT pose_id, pose_alias FROM pose 
            WHERE pose_alias IS NOT NULL
            AND pose_id IN {pset.str_ids}"""
            ).fetchall()

        else:
            records = self.execute(
                """SELECT pose_id, pose_alias FROM pose 
            WHERE pose_alias IS NOT NULL"""
            ).fetchall()

        d = {}
        for pose_id, pose_alias in records:
            d[pose_id] = pose_alias

        return d

    def get_pose_path_id_dict(self, pset: "PoseSet | None" = None) -> dict[str, int]:
        """Get a dictionary mapping :class:`.Pose` aliases to ID's"""

        if pset:
            records = self.execute(
                f"""
            SELECT pose_id, pose_path FROM pose 
            WHERE pose_path IS NOT NULL
            AND pose_id IN {pset.str_ids}"""
            ).fetchall()

        else:
            records = self.execute(
                """SELECT pose_id, pose_path FROM pose 
            WHERE pose_path IS NOT NULL"""
            ).fetchall()

        d = {}
        for pose_id, pose_path in records:
            d[pose_path] = pose_id

        return d

    def get_pose_id_obj_dict(self, pset: "PoseSet") -> "dict[id, Pose]":
        """Get a dictionary mapping :class:`.Pose` ID's to their objects"""

        query = "pose_id, pose_inchikey, pose_alias, pose_smiles, pose_reference, pose_path, pose_compound, pose_target, pose_mol, pose_fingerprint, pose_energy_score, pose_distance_score"
        records = self.select_where(
            query=query, table="pose", key=f"pose_id IN {pset.str_ids}", multiple=True
        )

        d = {}
        for entry in records:
            (
                pose_id,
                pose_inchikey,
                pose_alias,
                pose_smiles,
                pose_reference,
                pose_path,
                pose_compound,
                pose_target,
                pose_mol,
                pose_fingerprint,
                pose_energy_score,
                pose_distance_score,
            ) = entry
            d[pose_id] = Pose(
                self,
                pose_id,
                pose_inchikey,
                pose_alias,
                pose_smiles,
                pose_reference,
                pose_path,
                pose_compound,
                pose_target,
                pose_mol,
                pose_fingerprint,
                pose_energy_score,
                pose_distance_score,
            )
        return d

    def get_pose_id_interaction_tuples_dict(self, pset: "PoseSet") -> dict[int, set]:
        """Get a dictionary mapping :class:`.Pose` ID's to lists of `(interaction_type, feature_id)` tuples describing their interactions"""

        sql = f"""
        SELECT DISTINCT interaction_pose, feature_id, interaction_type FROM interaction 
        INNER JOIN feature ON interaction_feature = feature_id
        WHERE interaction_pose IN {pset.str_ids}
        """

        records = self.execute(sql).fetchall()

        ISETS = {}
        for pose_id, feature_id, interaction_type in records:
            values = ISETS.get(pose_id, set())
            values.add((interaction_type, feature_id))
            ISETS[pose_id] = values

        return ISETS

    def get_compound_id_inspiration_ids_dict(self) -> dict[int, set]:
        """Get a dictionary mapping :class:`.Compound` ID's to a set of :class:`Pose` ID's for the inspirations for the whole database"""

        sql = """
        SELECT compound_id, pose_id, inspiration_original FROM compound
        INNER JOIN scaffold ON compound_id = scaffold_base
        INNER JOIN pose ON compound_id = pose_compound
        INNER JOIN inspiration ON pose_id = inspiration_derivative
        """

        with mrich.spinner("Database.get_pose_id_interaction_ids_dict()"):
            records = self.execute(sql).fetchall()

        d = {}

        for compound_id, pose_id, inspiration_id in records:
            d[compound_id] = d.get(compound_id, set())
            d[compound_id].add(inspiration_id)

        return d

    def get_inspiration_tuples(self) -> list[int, int]:
        """Get a dictionary mapping :class:`.Pose` ID's to a set of :class:`Pose` ID's for the inspirations for the whole database"""
        sql = """SELECT inspiration_original, inspiration_derivative FROM inspiration"""
        return self.execute(sql).fetchall()

    def get_compound_id_obj_dict(self, cset: "CompoundSet") -> "dict[id, Compound]":
        """Get a dictionary mapping :class:`.Compound` ID's to their objects"""

        query = "compound_id, compound_inchikey, compound_alias, compound_smiles"
        records = self.select_where(
            query=query,
            table="compound",
            key=f"compound_id IN {cset.str_ids}",
            multiple=True,
        )

        d = {}
        for entry in records:
            compound_id, compound_inchikey, compound_alias, compound_smiles = entry
            d[compound_id] = Compound(
                self._animal,
                self,
                compound_id,
                compound_inchikey,
                compound_alias,
                compound_smiles,
                metadata=None,
                mol=None,
            )
        return d

    def get_interaction(self, *, id: int, table: str = "interaction") -> "Interaction":
        """Fetch the :class:`.Interaction` object with given ID

        :param id: the ID of the Interaction to retrieve
        :returns: :class:`.Interaction` object

        """

        from .interaction import Interaction

        result = self.select_all_where(table=table, key=f"interaction_id = {id}")

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
        ) = result

        return Interaction(
            db=self,
            id=id,
            feature_id=feature_id,
            pose_id=pose_id,
            type=type,
            family=family,
            atom_ids=atom_ids,
            prot_coord=prot_coord,
            lig_coord=lig_coord,
            distance=distance,
            angle=angle,
            energy=energy,
            table=table,
        )

    def get_reaction_map_from_products(
        self, product_ids: list[int]
    ) -> dict[tuple[str, int], set[int]]:
        """Get a dictionary mapping (reaction_type, product_id) tuples to sets of reactant_ids"""

        str_ids = str(tuple(product_ids)).replace(",)", ")")

        records = self.execute(
            f"""
            SELECT reaction_type, reaction_product, reaction_id, reactant_compound
            FROM reaction INNER JOIN reactant
            ON reaction_id = reactant_reaction
            WHERE reaction_product IN {str_ids}
        """
        ).fetchall()

        mapping = {}
        for reaction_type, reaction_product, reaction_id, reactant_compound in records:

            key = (reaction_type, reaction_product)

            if key not in mapping:
                mapping[key] = {}

            if reaction_id not in mapping[key]:
                mapping[key][reaction_id] = set()

            mapping[key][reaction_id].add(reactant_compound)

        return mapping

    def get_possible_reaction_ids(
        self,
        *,
        compound_ids: list[int],
    ) -> list[int]:
        """Given a set of reactant :class:`.Compound` ID's, compute which :class:`.Reaction` objects are possible (all reactants present).

        :param compound_ids: the list of reactant :class:`.Compound` IDs
        :returns: a list of :class:`.Reaction` ID's that are possible with the given reactants

        """

        compound_ids_str = str(tuple(compound_ids)).replace(",)", ")")

        result = self.execute(
            f"""
            WITH possible_reactants AS 
            (
                SELECT reactant_reaction, CASE WHEN reactant_compound IN {compound_ids_str} THEN reactant_compound END AS [possible_reactant] FROM reactant
            )

            , possible_reactions AS (
                SELECT reactant_reaction, COUNT(CASE WHEN possible_reactant IS NULL THEN 1 END) AS [count_null] FROM possible_reactants
                GROUP BY reactant_reaction
            )
            
            SELECT reactant_reaction FROM possible_reactions
            WHERE count_null = 0
        """
        ).fetchall()

        return [q for q, in result]

    def get_unsolved_reaction_tree(
        self,
        *,
        product_ids: list[int],
        debug: bool = False,
    ) -> "(CompoundSet, ReactionSet)":
        """Given a set of product :class:`.Compound` IDs, recursively solve for all the reactants (:class:`.CompoundSet`) and reactions (:class:`.ReactionSet`) that could be involved in their synthesis. N.B. This evaluates all synthesis branches.

        :param product_ids: list of product :class:`.Compound` IDs
        :param debug: increase verbosity for debugging, defaults to False
        :returns: a tuple of ``(reactants, reactions)``

        """

        from .cset import CompoundSet
        from .rset import ReactionSet

        all_reactants = set()
        all_reactions = set()

        # print(product_ids)

        # for product_id in product_ids:
        #     all_reactants.add(product_id)

        for i in range(300):

            if debug:
                mrich.var("recursive depth", i + 1)

            if debug:
                mrich.var("#products", len(product_ids))

            product_ids_str = str(tuple(product_ids)).replace(",)", ")")

            reaction_ids = self.select_where(
                table="reaction",
                query="DISTINCT reaction_id",
                key=f"reaction_product in {product_ids_str}",
                multiple=True,
                none="quiet",
            )

            reaction_ids = [q for q, in reaction_ids]

            if not reaction_ids:
                break

            for reaction_id in reaction_ids:
                all_reactions.add(reaction_id)

            if debug:
                mrich.var("#reactions", len(reaction_ids))

            reaction_ids_str = str(tuple(reaction_ids)).replace(",)", ")")

            reactant_ids = self.select_where(
                table="reactant",
                query="DISTINCT reactant_compound",
                key=f"reactant_reaction in {reaction_ids_str}",
                multiple=True,
            )

            if debug:
                mrich.var("#reactants", len(reactant_ids))

            reactant_ids = [q for q, in reactant_ids]

            if not reactant_ids:
                break

            for reactant_id in reactant_ids:
                all_reactants.add(reactant_id)

            product_ids = reactant_ids

        # all intermediates
        ids = self.execute(
            "SELECT DISTINCT reaction_product FROM reaction INNER JOIN reactant ON reaction.reaction_product = reactant.reactant_compound"
        ).fetchall()
        ids = [q for q, in ids]
        intermediates = CompoundSet(self, ids)

        # remove intermediates
        cset = CompoundSet(self, all_reactants)
        all_reactants = cset - intermediates

        # reactions
        all_reactions = ReactionSet(self, all_reactions)

        if debug:
            mrich.var("#all_reactants", len(all_reactants))
            mrich.var("#all_reactions", len(all_reactions))

        return all_reactants, all_reactions

    def get_reaction_price_estimate(
        self,
        *,
        reaction: Reaction,
    ) -> float:
        """Estimate the price of a :class:`.Reaction`

        :param reaction: :class:`.Reaction` object
        :returns: price estimate

        """

        # get reactants for a given reaction

        mrich.warning("Price estimate does not account for branching!")
        reactants, _ = self.get_unsolved_reaction_tree(
            product_ids=reaction.reactant_ids
        )

        # how to make sure that there are no branching reactions?!

        # sum lowest unit price for each reactant

        (price,) = self.execute(
            f"""
        WITH unit_prices AS 
        (
            SELECT quote_compound, MIN(quote_price/quote_amount) AS unit_price FROM quote 
            WHERE quote_compound IN {reactants.str_ids}
            GROUP BY quote_compound
        )
        SELECT SUM(unit_price) FROM unit_prices
        """
        ).fetchone()

        return price

    def get_possible_reaction_product_ids(
        self,
        *,
        reaction_ids: list[int],
    ) -> list[int]:
        """Given a set of :class:`.Reaction` IDs return the :class:`.Compound` IDs of their synthesis products

        :param reaction_ids: :class:`.Reaction` IDs
        :returns: list of :class:`.Compound` IDs

        """
        reaction_ids_str = str(tuple(reaction_ids)).replace(",)", ")")
        return [
            q
            for q, in self.select_where(
                query="DISTINCT reaction_product",
                table="reaction",
                key=f"reaction_id IN {reaction_ids_str}",
                multiple=True,
            )
        ]

    def get_subsite(self, *, id) -> "Subsite":
        """Get protein subsite with a given ID

        :param ID: the subsite ID
        :returns: :class:`.Subsite` object

        """

        from .subsite import Subsite

        name, target = self.select_where(
            table="subsite",
            key="id",
            value=id,
            multiple=False,
            query="subsite_name, subsite_target",
        )

        subsite = Subsite(db=self, id=id, name=name, target_id=target)

        return subsite

    def get_subsite_tag(self, *, id) -> "SubsiteTag":
        """Get subsite_tag with a given ID

        :param ID: the subsite_tag ID
        :returns: :class:`.SubsiteTag` object

        """

        from .subsite import SubsiteTag

        subsite_id, pose_id = self.select_where(
            table="subsite_tag",
            key="id",
            value=id,
            multiple=False,
            query="subsite_tag_ref, subsite_tag_pose",
        )

        subsite_tag = SubsiteTag(db=self, id=id, subsite_id=subsite_id, pose_id=pose_id)

        return subsite_tag

    def get_subsite_id(self, *, name: str, **kwargs) -> int | None:
        """Get protein Subsite ID with a given name

        :param name: the protein Subsite name
        :returns: the Subsite ID

        """

        table = "Subsite"
        entry = self.select_id_where(table=table, key="name", value=name, **kwargs)

        if entry:
            return entry[0]

        return None

    def get_subsite_name(self, *, id: str, **kwargs) -> int | None:
        """Get protein :class:`.Subsite` name with a given ID

        :param name: the protein :class:`.Subsite` ID
        :returns: the :class:`.Subsite` ID

        """

        table = "subsite"
        entry = self.select_where(
            query="subsite_name", table=table, key="id", value=id, **kwargs
        )

        if entry:
            return entry[0]

        return None

    def get_scaffold_similarity_dict(
        self, scaffolds: "CompoundSet | None" = None
    ) -> list[dict]:
        """Get a dictionary mapping scaffold :class:`.Compound` IDs to their superstructure's IDs"""

        sql = """
        SELECT scaffold_base as a, scaffold_superstructure as b, bfp_tanimoto(c.fp, d.fp) AS t
        FROM scaffold
        INNER JOIN compound_pattern_bfp AS c ON a = c.compound_id
        INNER JOIN compound_pattern_bfp AS d ON b = d.compound_id
        """

        if scaffolds:

            sql += f" WHERE a IN {scaffolds.str_ids}"

        records = self.execute(sql).fetchall()

        data = []
        for a, b, s in mrich.track(records):
            data.append(dict(base_id=a, superstructure_id=b, similarity=s))

        return data

    def get_reactant_product_tuples(
        self, compound_ids: list | None = None, deduplicated: bool = True
    ) -> set[tuple[int, int]]:
        """Get tuples of (reactant, product) :class:`.Compound` IDs"""

        sql = """
        SELECT reactant_compound, reaction_product FROM reactant
        INNER JOIN reaction ON reactant_reaction = reaction_id
        """

        if compound_ids:
            str_ids = str(tuple(compound_ids)).replace(",)", ")")
            sql += (
                f"WHERE reactant_compound IN {str_ids} OR reaction_product IN {str_ids}"
            )

        records = self.execute(sql)
        if deduplicated:
            return set((a, b) for a, b in records)
        else:
            return [(a, b) for a, b in records]

    def get_scaffold_tuples(
        self, compound_ids: list | None = None
    ) -> set[tuple[int, int]]:
        """Get tuples of (reactant, product) :class:`.Compound` IDs"""

        sql = """
        SELECT scaffold_base, scaffold_superstructure FROM scaffold
        """

        if compound_ids:
            str_ids = str(tuple(compound_ids)).replace(",)", ")")
            sql += f"WHERE scaffold_base IN {str_ids} OR scaffold_superstructure IN {str_ids}"

        records = self.execute(sql)
        return set((a, b) for a, b in records)

    ### COMPOUND QUERY

    def query_substructure(
        self,
        query: str,
        *,
        fast: bool = True,
        none: str = "error",
        smarts: bool = False,
    ) -> "CompoundSet":
        """Search for compounds by substructure

        :param query: SMILES string of the substructure
        :param fast: Use pattern binary fingerprint table to improve performance (Default value = True)
        :param none: define the behaviour for no matches, any value other than ``'error'`` will silently return empty data (Default value = 'error')
        :returns: :class:`.CompoundSet` object

        """

        if smarts:
            func = "mol_from_smarts"
        else:
            func = "mol_from_smiles"

        # smiles
        if isinstance(query, str):

            if fast:
                sql = f"""
                SELECT compound.compound_id, compound.compound_inchikey 
                FROM compound, compound_pattern_bfp AS bfp 
                WHERE compound.compound_id = bfp.compound_id 
                AND mol_is_substruct(compound.compound_mol, {func}(?))
                """

            else:
                sql = f"""
                SELECT compound_id, compound_inchikey FROM compound 
                WHERE mol_is_substruct(compound_mol, {func}(?))
                """

        else:

            raise NotImplementedError

        try:
            self.execute(sql, (query,))
        except sqlite3.OperationalError as e:
            mrich.var("sql", sql)
            raise

        result = self.cursor.fetchall()

        if not result and none == "error":
            mrich.error(f"No compounds with substructure {query}")
            return None

        from .cset import CompoundSet

        if not smarts:
            smiles = query
            try:
                smiles = sanitise_smiles(smiles, sanitisation_failed="error")
            except SanitisationError as e:
                mrich.error(f"Could not sanitise {smiles=}")
                mrich.error(str(e))
                return None
            except AssertionError:
                mrich.error(f"Could not sanitise {smiles=}")
                return None
                return c
            inchikey = inchikey_from_smiles(smiles)

            ids = [i for i, key in result if key != inchikey]
            return CompoundSet(self, ids)

        else:
            return CompoundSet(self, [i for i, _ in result])

    def query_similarity(
        self,
        query: str,
        threshold: float,
        return_similarity: bool = False,
        none="error",
    ) -> "CompoundSet | (CompoundSet, list[float])":
        """Search compounds by tanimoto similarity

        :param query: SMILES string
        :param threshold: similarity threshold to exceed
        :param return_similarity: return a list of similarity values together with the :class:`.CompoundSet` (Default value = False)
        :param none: define the behaviour for no matches, any value other than ``'error'`` will silently return empty data (Default value = 'error')
        :returns: :class:`.CompoundSet` and optionally a list of similarity values

        """

        from .cset import CompoundSet

        # smiles
        if isinstance(query, str):

            if return_similarity:
                sql = f"SELECT compound_id, bfp_tanimoto(mol_pattern_bfp(mol_from_smiles(?1), 2048), mol_pattern_bfp(compound.compound_mol, 2048)) as t FROM compound JOIN compound_pattern_bfp AS mfp USING(compound_id) WHERE mfp.compound_id match rdtree_tanimoto(mol_pattern_bfp(mol_from_smiles(?1), 2048), ?2) ORDER BY t DESC "
            else:
                sql = f"SELECT compound_id FROM compound_pattern_bfp AS bfp WHERE bfp.compound_id match rdtree_tanimoto(mol_pattern_bfp(mol_from_smiles(?1), 2048), ?2) "

        else:
            raise NotImplementedError

        try:
            self.execute(sql, (query, threshold))
        except sqlite3.OperationalError as e:
            mrich.var("sql", sql)
            raise

        result = self.cursor.fetchall()

        if not result and none == "error":
            mrich.error(f"No compounds with similarity >= {threshold} to {query}")
            return None

        if return_similarity:
            ids, similarities = zip(*result)
            cset = CompoundSet(self, ids)
            return cset, similarities

        ids = [r for r, in result]
        cset = CompoundSet(self, ids)

        return cset

    def query_exact(
        self,
        query: str,
        threshold: float = 0.989,
    ) -> "CompoundSet":
        """Search for exact match compounds (default similarity > 0.989)

        :param query: SMILES string
        :param threshold: similarity threshold to exceed

        """

        return self.query_similarity(query, 0.989, return_similarity=False)

    ### LOOKUPS / MAPPING

    def create_metadata_id_map(self, *, table: str, key: str) -> dict[str, int]:
        """Create a mapping between metadata[key] values to their respective parent record ID's

        :returns: dictionary mapping metadata[key] values to integer ID's

        """

        pairs = self.execute(
            f"""SELECT {table}_id, {table}_metadata FROM {table} WHERE {table}_metadata LIKE '%"{key}": "%'"""
        ).fetchall()
        from json import loads

        return dict(
            sorted(
                {loads(metadata)[key]: pose_id for pose_id, metadata in pairs}.items()
            )
        )

    ### COUNTING

    def count(
        self,
        table: str,
    ) -> int:
        """Count all entries in a table

        :param table: table to count entries from

        """

        sql = f"""SELECT COUNT(1) FROM {table}; """
        self.execute(sql)
        return self.cursor.fetchone()[0]

    def count_where(
        self,
        table: str,
        key: str,
        value=None,
    ):
        """Count all entries in a table where ``key==value``

        :param table: table to count entries from
        :param key: the key to match as ``{table}_{key} = {value}`` or the SQL string if ``value == None``
        :param value: the value to match (Default value = None)

        """
        if value is not None:
            where_str = f"{table}_{key} is {value}"
        else:
            where_str = key

        sql = f"""SELECT COUNT(1) FROM {table} WHERE {where_str};"""
        self.execute(sql)
        return self.cursor.fetchone()[0]

    ### ID SELECTION

    def min_id(
        self,
        table: str,
    ) -> int:
        """Get the minimal entry ID from a given table

        :param table: the database table to query
        :returns: the smallest entry ID

        """

        (id,) = self.select(table=table, query=f"MIN({table}_id)")
        return id

    def max_id(
        self,
        table: str,
    ) -> int:
        """Get the maximal entry ID from a given table

        :param table: the database table to query
        :returns: the largest entry ID

        """
        (id,) = self.select(table=table, query=f"MAX({table}_id)")
        return id

    def slice_ids(
        self,
        *,
        table: str,
        start: int | None,
        stop: int | None,
        step: int = 1,
        name: bool = False,
    ) -> list[int]:
        """Retrieve ID's matching a slice

        :param table: the database table to query
        :param start: return IDs equal to or larger than this value
        :param stop: return IDs smaller than this value
        :param step: return IDs in increments of this value (Default value = 1)
        :returns: matching IDs

        """

        min_id = self.min_id(table)
        max_id = self.max_id(table)

        start = start or min_id
        stop = stop or max_id + 1
        step = step or 1

        if not (start >= 0 and start <= max_id):
            raise IndexError(
                f"Slice {start=} outside of DB {table}_id range ({min_id}, {max_id})"
            )

        if not (stop >= 0 and stop <= max_id + 1):
            raise IndexError(
                f"Slice {stop=} outside of DB {table}_id range ({min_id}, {max_id})"
            )

        if step != 1:
            raise NotImplementedError(f"Slice {step=} not supported")

        ids = self.select_where(
            table=table,
            query=f"{table}_id",
            key=f"{table}_id >= {start} AND {table}_id < {stop}",
            multiple=True,
        )
        ids = [q for q, in ids]

        if name:
            return ids, f"{table}s[{start}:{stop}]"
        else:
            return ids

    ### PRUNING

    def prune_reactions(
        self,
        reactions: "ReactionSet",
    ) -> list[Reaction]:
        """Remove duplicate reactions

        :param reactions: :class:`.ReactionSet`
        :returns: list of pruned :class:`.Reaction` objects

        """

        pruned = []
        del_list = []

        for i, reaction in enumerate(reactions):

            matches = [r for r in pruned if r == reaction]

            if not matches:
                pruned.append(reaction)

            else:
                del_list.append(reaction)

        for reaction in del_list:
            mrich.warning(f"Deleted duplicate {reaction=}")
            self.delete_where("reaction", "id", reaction.id)

        return pruned

    def remove_metadata_list_item(
        self,
        *,
        table: str,
        key: str,
        value,
        remove_empty: bool = True,
    ) -> None:
        """Remove a specific item from list-like values associated with a given key from all metadata entries in a given table

        :param table: the database table to query
        :param key: the :class:`.Metadata` key to match
        :param value: the value to remove from the list
        :param remove_empty: remove the key from the metadata if the list is empty (Default value = True)

        """

        # get id's with specific metadata key and value
        value_str = json.dumps(value)
        result = self.select_where(
            query=f"{table}_id, {table}_metadata",
            table=table,
            key=f"{table}_metadata LIKE '%\"export\": [%{value_str}%]%'",
            multiple=True,
            none="quiet",
        )

        # loop over all matches
        for id, metadata_str in result:

            # read the metadata
            metadata = json.loads(metadata_str)

            # modify the metadata
            index = metadata[key].index(value)
            metadata[key].pop(index)

            if remove_empty and not metadata[key]:
                del metadata[key]

            # update the database
            metadata_str = json.dumps(metadata)
            self.update(
                table=table,
                id=id,
                key=f"{table}_metadata",
                value=metadata_str,
                commit=False,
            )

        # only runs if non-zero matches
        else:
            # commit the changes
            self.commit()

    ### TABLE INFO

    def print_table(
        self,
        table: str,
    ) -> None:
        """Print a table's entries

        :param table: the table to print

        """

        # mrich.print(self.cursor.fetchall())

        from rich.table import Table

        tab = Table()

        for col in self.column_names(table):
            tab.add_column(col.removeprefix(table).removeprefix("_"))

        self.execute(f"SELECT * FROM {table}")
        for record in self.cursor.fetchall():
            record = [str(v) for v in record]
            tab.add_row(*record)

        mrich.print(tab)

    def table_info(
        self,
        table: str,
    ) -> list[tuple]:
        """Print a table's schema

        :param table: the table to print

        """

        self.execute(f"PRAGMA table_info({table})")
        return self.cursor.fetchall()

    def column_names(self, table: str) -> list[str]:
        """Get the column names of the given table"""
        table_info = self.table_info(table)
        return [i[1] for i in table_info]

    ### DUNDERS

    def __str__(self):
        """Unformatted string representation"""
        return f"Database @ {self.path.resolve()}"

    def __repr__(self):
        """ANSI Formatted string representation"""
        return f"{mcol.bold}{mcol.underline}{self}{mcol.clear}"

    def __rich__(self) -> str:
        """Representation for mrich"""
        return f"[bold underline]{self}"


class LegacyDatabaseError(Exception): ...
