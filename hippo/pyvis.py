import mrich


def get_scaffold_network(
    animal,
    compounds="CompoundSet | None",
    scaffolds="CompoundSet | None",
    # filename: "str | Path" = "network.html",
    notebook: bool = True,
    depth: int = 5,
    scaffold_tag: str | None = None,
    exclude_tag: str | None = None,
    physics: bool = True,
    arrows: bool = True,
) -> "pyvis.network.Network":
    """Use PyVis to display a network of molecules connected by scaffold relationships in the database"""

    from pyvis.network import Network
    from molparse.rdkit import smiles_to_pngstr

    net = Network(notebook=notebook, cdn_resources="in_line")

    nodes = set()
    edges = set()

    arrows = "to" if arrows else None

    def add_node(compound):

        if compound.id in nodes:
            return

        pngstr = smiles_to_pngstr(compound.smiles)

        net.add_node(
            compound.id,
            label=compound.alias or str(compound),
            title=str(compound),
            shape="circularImage",
            image=f"data:image/png;base64,{pngstr}",
            physics=physics,
        )

        nodes.add(compound.id)

    def add_edge(base, compound):

        key = (base.id, compound.id)

        if key in edges:
            return

        net.add_edge(
            base.id,
            compound.id,
            arrows=arrows,
            # label="PDE5",
            # title="Protein target",
            # color="purple"
        )

        edges.add(key)

    def get_scaffold_records(scaffolds=None, compounds=None):
        if scaffolds:
            return animal.db.select_all_where(
                table="scaffold",
                key=f"scaffold_base IN {scaffolds.str_ids}",
                multiple=True,
                none="quiet",
            )
        elif compounds:
            return animal.db.select_all_where(
                table="scaffold",
                key=f"scaffold_superstructure IN {compounds.str_ids}",
                multiple=True,
                none="quiet",
            )
        raise ValueError

    if compounds and not scaffolds:

        mrich.var("recursion depth", depth)
        mrich.var("#compounds", len(compounds))

        records = []
        n = 0
        while n < depth:

            if not compounds:
                break

            results = get_scaffold_records(compounds=compounds)

            if results:
                records.extend(results)
                compounds = animal.compounds[[a for a, b in results]]
                n += 1
            else:
                break

        if n == depth:
            mrich.warning(
                "Reached recursion depth. More scaffolds may be in the database"
            )

    elif scaffolds and not compounds:

        mrich.var("recursion depth", depth)
        mrich.var("#scaffolds", len(scaffolds))

        records = []
        n = 0
        while n < depth:

            if not scaffolds:
                break

            results = get_scaffold_records(scaffolds=scaffolds)

            if results:
                records.extend(results)
                scaffolds = animal.compounds[[b for a, b in results]]
                n += 1
            else:
                break

        if n == depth:
            mrich.warning(
                "Reached recursion depth. More superstructures may be in the database"
            )

    else:
        raise ValueError

    mrich.var("#edges", len(records))

    if records:
        for base_id, compound_id in mrich.track(
            records, prefix="Adding nodes and edges"
        ):

            base = animal.db.get_compound(id=base_id)

            base_tags = base.tags

            if scaffold_tag and scaffold_tag not in base_tags:
                continue

            compound = animal.db.get_compound(id=compound_id)

            if exclude_tag and (
                exclude_tag in base_tags or exclude_tag in compound.tags
            ):
                continue

            add_node(base)
            add_node(compound)
            add_edge(base, compound)

    return net
