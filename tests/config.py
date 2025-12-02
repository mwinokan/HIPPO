## CONFIGURE TESTING DATA
TARGET = "SARS2_Nprot"
PROPOSAL = "lb32627-93"
STACK = "production"

DB = dict(
    username="postgres",
    password="hippo",
    host="localhost",
    port=5432,
)

TARGET = "SARS2_Nprot"

## DISABLE TESTS

CLEANUP = True
DOWNLOAD = False
SETUP = True
ADD_HITS = True
SCAFFOLDS = True
SUBSITES = True

## CONFIGURE CLEANUP

CLEANUP_FILES = [
    # DB,
    # f"{TARGET}.tar.gz",
]

CLEANUP_DIRS = [
    # TARGET,
]
