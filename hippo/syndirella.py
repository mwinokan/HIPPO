def reactions_from_row(*, animal, row, num_steps):

    reaction_ids = set()

    for step in range(num_steps):

        step += 1

        # get relevant fields
        reaction_name = row[f"{step}_reaction"]

        product = None
        reactants = []

        if reactant1_smiles := row[f"{step}_r1_smiles"]:
            reactants.append(animal.register_compound(smiles=reactant1_smiles))

        if reactant2_smiles := row[f"{step}_r2_smiles"]:
            reactants.append(animal.register_compound(smiles=reactant2_smiles))

        if product_smiles := row[f"{step}_product_smiles"]:
            product = animal.register_compound(smiles=product_smiles)

        assert product
        assert reactants

        reaction = animal.register_reaction(
            type=reaction_name,
            product=product,
            reactants=reactants,
        )

        assert reaction

        reaction_ids.add(reaction.id)

    return animal.reactions[reaction_ids]
