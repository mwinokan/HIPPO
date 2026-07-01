## CONFIGURE TESTING DATA
TARGET = 'SARS2_Nprot'
PROPOSAL = 'lb32627-93'
STACK = 'production'

## CONFIGURE CLEANUP

CLEANUP_FILES = [
    f'{TARGET}.tar.gz',
]

CLEANUP_DIRS = [
    TARGET,
]

## CONFIGURE DATABASE

## CONFIGURE TESTS

CLEANUP = True
DOWNLOAD = True
SETUP = True
ADD_HITS = True
SCAFFOLDS = True
SUBSITES = True

### SQLITE

DB = 'db_test.sqlite'

CLEANUP_FILES.append(DB)

### POSTGRES

# SCAFFOLDS = False

# local testing

# from os import environ

# DB = dict(
#     username=environ["HIPPO_POSTGRES_USERNAME"],
#     password=environ["HIPPO_POSTGRES_PASSWORD"],
#     host="localhost",
#     port=5555,
#     dbname="postgres",
# )
