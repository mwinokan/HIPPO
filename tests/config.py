## CONFIGURE TESTING DATA
TARGET = "SARS2_Nprot"
PROPOSAL = "lb32627-93"
STACK = "production"

## CONFIGURE CLEANUP

CLEANUP_FILES = [
    # f"{TARGET}.tar.gz",
]

CLEANUP_DIRS = [
    # TARGET,
]

## CONFIGURE DATABASE

### SQLITE

DB = "db_test.sqlite"

CLEANUP_FILES.append(DB)

### POSTGRES

DB = dict(
    username="postgres",
    password="hippo",
    host="localhost",
    port=5432,
)

## DISABLE TESTS

CLEANUP = True
DOWNLOAD = False
SETUP = True
ADD_HITS = True
SCAFFOLDS = True
SUBSITES = True
