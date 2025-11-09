from .animal import HIPPO
import mrich
from typer import Typer
from pathlib import Path

app = Typer()


def setup_animal(
    database: str,
    backup: bool = True,
    update_legacy: bool = False,
) -> HIPPO:
    animal = HIPPO("CLI", database, update_legacy=update_legacy)
    if backup:
        animal.db.backup()
    return animal


@app.command()
def backup(database: str):
    """Backup database file"""
    mrich.h1("hippo.backup")

    mrich.h3("Params")
    mrich.var("database", database)

    from hippo.db import backup

    backup(database)
    mrich.success("Successfully backed up")


@app.command()
def update_legacy(database: str, backup: bool = True):
    """Update legacy database format"""

    mrich.h1("hippo.update_legacy")

    mrich.h3("Params")
    mrich.var("database", database)

    if backup:
        from hippo.db import backup

        backup(database)

    animal = setup_animal(database, backup=False, update_legacy=True)
    mrich.success("Successfully updated database format")


@app.command()
def calculate_scaffolds(
    database: str,
    backup: bool = True,
):
    """Calculate scaffold/superstructure relationships for all compounds"""

    mrich.h1("hippo.calculate_scaffolds")

    mrich.h3("Params")
    mrich.var("database", database)
    mrich.var("backup", backup)

    mrich.h3("Animal")
    animal = setup_animal(database=database, backup=backup)

    mrich.h3("State Before")
    mrich.var("scaffolds", animal.scaffolds)
    mrich.var("elabs", animal.elabs)

    mrich.h3("Calculation")
    animal.db.calculate_all_scaffolds()

    mrich.h3("State After")
    mrich.var("scaffolds", animal.scaffolds)
    mrich.var("elabs", animal.elabs)

    mrich.success("Completed")


@app.command()
def calculate_interactions(
    database: str,
    prolif: bool = False,
    backup: bool = True,
    force: bool = False,
):
    """Calculate interactions for all poses"""

    mrich.h1("hippo.calculate_interactions")

    mrich.h3("Params")
    mrich.var("database", database)
    mrich.var("backup", backup)

    mrich.h3("Animal")
    animal = setup_animal(database=database, backup=backup)

    mrich.h3("State Before")
    mrich.var("#fingerprinted", animal.poses.num_fingerprinted)

    mrich.h3("Calculation")

    n = len(animal.poses)
    for i, pose in mrich.track(enumerate(animal.poses), total=n):

        mrich.set_progress_prefix(f"{i}/{n}")

        try:
            if prolif:
                pose.calculate_prolif_interactions(force=force)
            else:
                pose.calculate_interactions(force=force)

        except Exception as e:
            mrich.error(e)
            mrich.error("Could not fingerprint pose")
            continue

    mrich.h3("State After")
    mrich.var("#fingerprinted", animal.poses.num_fingerprinted)

    mrich.success("Completed")


@app.command()
def verify():
    """Verify installation"""

    import os

    file_path = "_test.sqlite"

    try:
        animal = setup_animal(file_path, backup=False)
        c = animal.register_compound(smiles="COc1ccc2sc(N)nc2c1")
        c.mol
        mrich.success("HIPPO/rdkit/chemicalite installations are compatible")
    except Exception as e:
        mrich.error(e)

    if os.path.exists(file_path):
        os.remove(file_path)


@app.command()
def tag_summary(
    database: str,
):
    """Print a table of statistics for all tags in the database"""

    mrich.h1("hippo.tag_summary")

    mrich.h3("Params")
    mrich.var("database", database)

    animal = setup_animal(database=database, backup=False)
    animal.tags.summary()


def main():
    app()


if __name__ == "__main__":
    main()
