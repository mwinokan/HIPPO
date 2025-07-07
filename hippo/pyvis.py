import mrich


def get_scaffold_network(
    animal,
    compounds="CompoundSet | None",
    scaffolds="CompoundSet | None",
    # filename: "str | Path" = "network.html",
    notebook: bool = True,
    depth: int = 1,
    scaffold_tag: str | None = None,
) -> "pyvis.network.Network":
    """Use PyVis to display a network of molecules connected by scaffold relationships in the database"""

    from pyvis.network import Network
    from molparse.rdkit import smiles_to_pngstr

    net = Network(notebook=notebook)

    if depth != 1:
        raise NotImplementedError

    nodes = set()
    edges = set()

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
        )

        nodes.add(compound.id)

    def add_edge(base, compound):

        key = (base.id, compound.id)

        if key in edges:
            return

        net.add_edge(
            base.id,
            compound.id,
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
            )
        elif compounds:
            return animal.db.select_all_where(
                table="scaffold",
                key=f"scaffold_superstructure IN {compounds.str_ids}",
                multiple=True,
            )
        raise ValueError

    if compounds and not scaffolds:
        mrich.var("#superstructures", len(compounds))
        records = get_scaffold_records(compounds=compounds)

    elif scaffolds and not compounds:
        mrich.var("#scaffolds", len(scaffolds))
        records = get_scaffold_records(scaffolds=scaffolds)

    else:
        raise ValueError

    mrich.var("#edges", len(records))

    if records:
        for base_id, compound_id in mrich.track(
            records, prefix="Adding scaffolds and edges"
        ):

            base = animal.db.get_compound(id=base_id)

            if scaffold_tag and scaffold_tag not in base.tags:
                continue

            compound = animal.db.get_compound(id=compound_id)

            add_node(base)
            add_node(compound)
            add_edge(base, compound)

    return net
