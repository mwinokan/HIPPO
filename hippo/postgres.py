"""PostgreSQL database wrapper class using psycopg3"""

import mcol
import mrich

import psycopg

from .db import Database


class PostgresDatabase(Database):
    """Wrapper to connect to a HIPPO Postgres database.

    .. attention::

            :class:`.PostGresDatabase` objects should not be created directly. Instead use the methods in :class:`.HIPPO` to interact with data in the database. See :doc:`getting_started` and :doc:`insert_elaborations`.

    """

    SQL_STRING_PLACEHOLDER = "%s"
    SQL_PK_DATATYPE = "SERIAL"
    SQL_SCHEMA = "hippo"
    SQL_SCHEMA_PREFIX = f"{SQL_SCHEMA}."

    ERROR_UNIQUE_VIOLATION = psycopg.errors.UniqueViolation

    SQL_CREATE_TABLE_COMPOUND = """CREATE TABLE hippo.compound(
        compound_id SERIAL PRIMARY KEY,
        compound_inchikey TEXT,
        compound_alias TEXT,
        compound_smiles TEXT,
        compound_base INTEGER,
        compound_mol MOL,
        compound_pattern_bfp bit(2048),
        compound_morgan_bfp bit(2048),
        compound_metadata TEXT,
        FOREIGN KEY (compound_base) REFERENCES hippo.compound(compound_id),
        CONSTRAINT UC_compound_inchikey UNIQUE (compound_inchikey),
        CONSTRAINT UC_compound_alias UNIQUE (compound_alias),
        CONSTRAINT UC_compound_smiles UNIQUE (compound_smiles)
    );
    """

    SQL_CREATE_TABLE_POSE = """CREATE TABLE hippo.pose(
        pose_id SERIAL PRIMARY KEY,
        pose_inchikey TEXT,
        pose_alias TEXT,
        pose_smiles TEXT,
        pose_reference INTEGER,
        pose_path TEXT,
        pose_compound INTEGER,
        pose_target INTEGER,
        pose_mol MOL,
        pose_fingerprint INTEGER,
        pose_energy_score REAL,
        pose_distance_score REAL,
        pose_inspiration_score REAL,
        pose_metadata TEXT,
        FOREIGN KEY (pose_compound) REFERENCES hippo.compound(compound_id),
        CONSTRAINT UC_pose_alias UNIQUE (pose_alias),
        CONSTRAINT UC_pose_path UNIQUE (pose_path)
    );
    """

    SQL_INSERT_COMPOUND = """
    INSERT INTO hippo.compound(
        compound_inchikey, 
        compound_smiles, 
        compound_mol, 
        -- compound_pattern_bfp, 
        -- compound_morgan_bfp, 
        compound_alias
    )
    VALUES(
        %(inchikey)s, 
        %(smiles)s, 
        mol_from_smiles(%(smiles)s), 
        -- mol_pattern_bfp(mol_from_smiles(%(smiles)s), 2048), 
        -- mol_morgan_bfp(mol_from_smiles(%(smiles)s), 2, 2048), 
        %(alias)s
    )
    RETURNING compound_id;
    """

    def __init__(
        self,
        animal: "HIPPO",
        username: str,
        password: str,
        host: str = "localhost",
        port: int = 5432,
        dbname: str = "hippo",
        update_legacy: bool = False,
        auto_compute_bfps: bool = False,
        create_blank: bool = True,
        check_legacy: bool = False,
        create_indexes: bool = True,
        update_indexes: bool = False,
        debug: bool = True,
    ) -> None:
        """PostgresDatabase initialisation"""

        assert isinstance(username, str)
        assert isinstance(password, str)
        assert isinstance(port, int)

        if debug:
            mrich.debug("hippo.PostgresDatabase.__init__()")

        self._username = username
        self._password = password
        self._port = port
        self._host = host

        self._connection = None
        self._cursor = None
        self._animal = animal
        self._auto_compute_bfps = auto_compute_bfps
        self._engine = "psycopg"
        self._dbname = dbname

        if debug:
            mrich.debug(f"PostgresDatabase.username = {self.username}")
            mrich.debug(f"PostgresDatabase.password = {self.password}")
            mrich.debug(f"PostgresDatabase.host = {self.host}")
            mrich.debug(f"PostgresDatabase.port = {self.port}")

        self.connect()

        if not self.table_names:

            if create_blank:
                self.execute("CREATE SCHEMA IF NOT EXISTS hippo;")
                self.create_blank_db()
            else:
                mrich.error("Database is empty!", self.path)
                raise ValueError(
                    "Database is empty! Check connection or run with create_blank=True"
                )

        if check_legacy:
            self.check_schema(update=update_legacy)

        if create_indexes:
            self.create_indexes(update=update_indexes, debug=debug)

    ### PROPERTIES

    @property
    def path(self) -> None:
        """PostgresDatabase path"""
        # raise NotImplementedError("PostgresDatabase has no path")
        return f"postgresql://{self.username}:{self.password}@{self.host}:{self.port}"

    @property
    def username(self) -> str:
        """PostgresDatabase username"""
        return self._username

    @property
    def dbname(self) -> str:
        """PostgresDatabase dbname"""
        return self._dbname

    @property
    def password(self) -> str:
        """PostgresDatabase password"""
        return self._password

    @property
    def host(self) -> str:
        """PostgresDatabase host"""
        return self._host

    @property
    def port(self) -> int:
        """PostgresDatabase port"""
        return self._port

    @property
    def table_names(self) -> list[str]:
        """List of all the table names in the database"""
        results = self.execute(
            f"""
            SELECT table_name
            FROM information_schema.tables
            WHERE table_schema = '{self.SQL_SCHEMA}'
            AND table_type = 'BASE TABLE';
        """
        ).fetchall()
        return [n for n, in results]

    def index_names(self) -> list[str]:
        """Get the index names"""

        cursor = self.execute(
            """
            SELECT indexname
            FROM pg_indexes
            WHERE schemaname = 'public';
        """
        )

        return [n for n, in cursor]

    @property
    def total_changes(self) -> int:
        """Return the current transaction ID as a proxy of sqlite's total_changes."""
        cursor = self.execute("SELECT txid_current()")
        return cursor.fetchone()[0]

    ### GENERAL SQL

    def connect(self, debug: bool = True) -> None:
        """Connect to the database"""

        if debug:
            mrich.debug("hippo.PostgresDatabase.connect()")

        conn = None

        try:
            conn = psycopg.connect(
                user=self.username,
                host=self.host,
                password=self.password,
                port=self.port,
                dbname=self.dbname,
            )

        except Exception as e:
            mrich.error("Could not connect to", self.path)
            mrich.error(e)
            raise

        self._connection = conn
        self._cursor = conn.cursor()

    def execute(
        self, sql, payload=None, *, retry: float | None = 1, debug: bool = False
    ):
        """Execute arbitrary SQL with retry if database is locked."""
        if debug:
            mrich.debug(sql)

        # while True:
        try:
            if payload:
                return self.cursor.execute(sql, payload)
            else:
                return self.cursor.execute(sql)
        # except sqlite3.OperationalError as e:
        #     if "database is locked" in str(e) and retry:
        #         with mrich.clock(
        #             f"SQLite Database is locked, waiting {retry} second(s)..."
        #         ):
        #             time.sleep(retry)
        #         mrich.print("[debug]SQLite Database was locked, retrying...")
        #         continue  # retry without recursion
        #     elif "syntax error" in str(e):
        #         mrich.error(sql)
        #         mrich.error(payload)
        #         raise
        #     else:
        #         raise
        except Exception as e:
            # mrich.print(sql)
            # mrich.print(payload)
            raise

    def rollback(self) -> None:
        """rollback (not relevant for sqlite)"""
        self.connection.rollback()

    def sql_return_id_str(self, key: str) -> str:
        """Add this to SQL queries to return the entry primary key"""
        return f"RETURNING {key}_id"

    def get_lastrowid(self) -> int:
        """Get ID of last inserted row"""
        return self.cursor.fetchone()[0]

    ### CREATE TABLES

    def create_table_pattern_bfp(self) -> None:
        """Create the pattern_bfp table"""
        mrich.warning(
            "HIPPO.PostgresDatabase.create_table_pattern_bfp(): NotImplemented"
        )

        return

        mrich.debug("HIPPO.PostgresDatabase.create_table_pattern_bfp()")

        sql = """
        CREATE VIRTUAL TABLE compound_pattern_bfp 
        USING rdtree(compound_id, fp bits(2048))
        """

        self.execute(sql)

    ### METHODS

    def _clear_schema(self) -> None:
        """Empty the Database schema entirely and recreate it"""

        self.execute(
            f"""
            DROP SCHEMA IF EXISTS {self.SQL_SCHEMA} CASCADE;
            CREATE SCHEMA {self.SQL_SCHEMA};
        """
        )

        self.commit()

    ### DUNDERS

    def __str__(self):
        """Unformatted string representation"""
        return f"Database @ {self.path}"
