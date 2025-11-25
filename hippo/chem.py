"""functions for validating chemistry"""

import mrich

"""

Checks
======

- Num heavy atoms difference
- Formula checks
- Num rings difference

"""

SUPPORTED_CHEMISTRY = {
    "Amidation": {
        "heavy_atoms_diff": 1,
        "rings_diff": 0,
        "atomtype": {
            "removed": {"O": 1, "H": 2},
        },
    },
    "Ester_amidation": {
        "heavy_atoms_diff": ">=3",
        "rings_diff": 0,
        "atomtype": {
            "removed": {"O": ">=1", "*": "*"},
        },
    },
    "Williamson_ether_synthesis": {
        "heavy_atoms_diff": 1,
        "rings_diff": 0,
        "atomtype": {
            "removed": {"Ha": 1, "H": 1},  # any halogen
        },
    },
    "N-Boc_deprotection": {
        "heavy_atoms_diff": 7,
        "rings_diff": 0,
        "atomtype": {
            "removed": {"O": 2, "C": 5, "H": 8},
        },
    },
    "TBS_alcohol_deprotection": {
        "heavy_atoms_diff": 7,
        "rings_diff": 0,
        "atomtype": {
            "removed": {"C": 6, "Si": 1, "H": 14},
        },
    },
    "Sp3-sp2_Suzuki_coupling": {
        # "heavy_atoms_diff": 10,
        "heavy_atoms_diff": ">=4",
        "rings_diff": ">=0",
        "atomtype": {
            # "removed": {"C": 6, "O": 2, "B": 1, "Ha": 1, "H": 12},  # any halogen
            "removed": {"C": ">=0", "O": 2, "B": 1, "Ha": 1, "H": ">=2"},  # any halogen
        },
    },
    "Sp2-sp2_Suzuki_coupling": {
        "heavy_atoms_diff": ">=4",
        "rings_diff": ">=0",
        "atomtype": {
            "removed": {"C": ">=0", "O": 2, "B": 1, "Ha": 1, "H": ">=2"},  # any halogen
        },
    },
    "Buchwald-Hartwig_amidation_with_amide-like_nucleophile": {
        "heavy_atoms_diff": 1,
        "rings_diff": 0,
        "atomtype": {
            "removed": {"Ha": 1, "H": 1},
        },
    },
    "Buchwald-Hartwig_amination": {
        "heavy_atoms_diff": 1,
        "rings_diff": 0,
        "atomtype": {
            "removed": {"Ha": 1, "H": 1},
        },
    },
    "Nucleophilic_substitution_with_amine": {
        "heavy_atoms_diff": 1,
        "rings_diff": 0,
        "atomtype": {
            "removed": {"Ha": 1, "H": 1},
        },
    },
    "N-nucleophilic_aromatic_substitution": {
        "heavy_atoms_diff": 1,
        "rings_diff": 0,
        "atomtype": {
            "removed": {"Ha": 1, "H": 1},  # any halogen
        },
    },
    "Reductive_amination": {
        "heavy_atoms_diff": 1,
        "rings_diff": 0,
        "atomtype": {
            "removed": {"O": 1},  # any halogen
        },
    },
    "Mitsunobu_reaction_with_amine_alcohol_and_thioalcohol": {
        "heavy_atoms_diff": 1,
        "rings_diff": 0,
        "atomtype": {
            "removed": {"O": 1, "H": ">=1"},
        },
    },
    "Steglich_esterification": {
        "heavy_atoms_diff": 1,
        "rings_diff": 0,
        "atomtype": {
            "removed": {"O": 1, "H": 2},
        },
    },
    "Benzyl_alcohol_deprotection": {
        "heavy_atoms_diff": 7,
        "rings_diff": 1,
        "atomtype": {
            "removed": {"C": 7, "H": 6},
        },
    },
    "N-Bn_deprotection": {
        "heavy_atoms_diff": 7,
        "rings_diff": 1,
    },
    "Formation_of_urea_from_two_amines": {
        "heavy_atoms_diff": -2,
        "rings_diff": 0,
    },
    "Amide_Schotten-Baumann_with_amine": {
        "heavy_atoms_diff": 1,
        "rings_diff": 0,
    },
    "Nucleophilic_substitution": {
        "heavy_atoms_diff": 1,
        "rings_diff": 0,
    },
}


def check_reaction_types(types: list[str]) -> None:
    """
    Prints a warning if any of the reaction type strings in ``types`` are not in ``SUPPORTED_CHEMISTRY``

    :param types: A list of reaction type strings to check
    """

    for reaction_type in types:
        if reaction_type not in SUPPORTED_CHEMISTRY:
            mrich.error(f"Can't check chemistry of unsupported {reaction_type=}")


def check_chemistry(
    reaction_type: str,
    reactants: "CompoundSet",
    product: "Compound",
    debug: bool = False,
) -> bool:
    """Check chemistry of given reaction"""

    if reaction_type not in SUPPORTED_CHEMISTRY:

        mrich.var("reactants", reactants.ids)
        mrich.var("product", product)

        raise UnsupportedChemistryError(f"Unsupported {reaction_type=}")

    assert reactants
    assert product

    CHEMISTRY = SUPPORTED_CHEMISTRY[reaction_type]

    if "heavy_atoms_diff" in CHEMISTRY:
        check = check_count_diff(
            "heavy_atoms", reaction_type, reactants, product, debug=debug
        )
        if not check:
            return False

    if "rings_diff" in CHEMISTRY:
        check = check_count_diff(
            "rings", reaction_type, reactants, product, debug=debug
        )
        if not check:
            return False

    if "atomtype" in CHEMISTRY:
        check = check_atomtype_diff(reaction_type, reactants, product, debug=debug)
        if not check:
            return False

    if debug:
        mrich.success(f"{reaction_type}: All OK")

    return True


def check_count_diff(
    check_type: str,
    reaction_type: str,
    reactants: "CompoundSet",
    product: "Compound",
    debug: bool = False,
):
    """Check integer difference"""

    # get target value
    diff = SUPPORTED_CHEMISTRY[reaction_type][f"{check_type}_diff"]

    # get attribute name
    attr = f"num_{check_type}"

    # get values
    reac_count = getattr(reactants, attr)
    prod_count = getattr(product, attr)
    if debug:
        mrich.var(f"#{check_type} reactants", reac_count)
    if debug:
        mrich.var(f"#{check_type} product", prod_count)

    # check against target value
    if isinstance(diff, str):

        assert diff.startswith(">="), diff

        diff = int(diff[2:])

        if reac_count - prod_count < diff:
            if debug:
                mrich.error(
                    f"{reaction_type}: #{check_type} {(reac_count - prod_count)=} FAIL"
                )
            return False

        elif debug:
            mrich.success(f"{reaction_type}: #{check_type} OK")

    else:

        if reac_count - diff != prod_count:
            if debug:
                mrich.error(f"{reaction_type}: #{check_type} FAIL")
            return False

        elif debug:
            mrich.success(f"{reaction_type}: #{check_type} OK")

    return True


def check_atomtype_diff(
    reaction_type: str,
    reactants: "CompoundSet",
    product: "Compound",
    debug: bool = False,
) -> bool:
    """check atomtypes"""

    check_type = "atomtype"

    # get values
    reac = reactants.atomtype_dict
    prod = product.atomtype_dict

    if debug:
        mrich.var("reactants.atomtype_dict", str(reac))
        mrich.var("product.atomtype_dict", str(prod))

    if "removed" in SUPPORTED_CHEMISTRY[reaction_type]["atomtype"]:
        removal = check_specific_atomtype_diff(
            reaction_type, prod, reac, removal=True, debug=debug
        )

        if not removal:
            return False

    if "added" in SUPPORTED_CHEMISTRY[reaction_type]["atomtype"]:
        addition = check_specific_atomtype_diff(
            reaction_type, prod, reac, removal=False, debug=debug
        )

        if not addition:
            return False

    if debug:
        mrich.success(f"{reaction_type}: atomtypes OK")

    return True


def check_specific_atomtype_diff(
    reaction_type: str,
    prod: "Compound",
    reac: "Compound",
    removal: bool = False,
    debug: bool = False,
) -> bool:
    """check specific atomtype difference"""

    if removal:
        add_str = "removed"
    else:
        add_str = "added"

    add_dict = SUPPORTED_CHEMISTRY[reaction_type]["atomtype"][add_str]

    if not add_dict:
        return True

    if debug:
        mrich.var(add_str, str(add_dict))

    for symbol, count in add_dict.items():

        if symbol == "Ha":
            p_count = halogen_count(prod)
            r_count = halogen_count(reac)

        elif symbol == "*":
            assert count == "*", (symbol, count)
            if debug:
                mrich.debug("Allowing wildcard atomtype differences")
            continue

        else:
            p_count = prod[symbol] if symbol in prod else 0
            r_count = reac[symbol] if symbol in reac else 0

        if isinstance(count, str):

            assert count.startswith(">="), (symbol, count)

            count = int(count[2:])

            if removal and r_count - p_count < count:
                if debug:
                    mrich.error(
                        f"{symbol}: {r_count=} - {p_count=} >= {r_count - p_count}"
                    )
                    mrich.error(
                        f"{reaction_type}: atomtype removal {symbol} × {count} FAIL"
                    )
                return False

            elif not removal and p_count - r_count < count:
                if debug:
                    mrich.error(
                        f"{symbol}: {p_count=} - {r_count=} >= {p_count - r_count}"
                    )
                    mrich.error(
                        f"{reaction_type}: atomtype addition {symbol} × {count} FAIL"
                    )
                return False

        else:

            if removal and r_count - p_count != count:
                if debug:
                    mrich.error(
                        f"{symbol}: {r_count=} - {p_count=} = {r_count - p_count}"
                    )
                    mrich.error(
                        f"{reaction_type}: atomtype removal {symbol} × {count} FAIL"
                    )
                return False

            elif not removal and p_count - r_count != count:
                if debug:
                    mrich.error(
                        f"{symbol}: {p_count=} - {r_count=} = {p_count - r_count}"
                    )
                    mrich.error(
                        f"{reaction_type}: atomtype addition {symbol} × {count} FAIL"
                    )
                return False

    return True


def halogen_count(atomtype_dict: dict[str, int]) -> int:
    """Count halogens"""
    count = 0
    symbols = ["F", "Cl", "Br", "I"]
    for symbol in symbols:
        if symbol in atomtype_dict:
            count += atomtype_dict[symbol]
    return count


class InvalidChemistryError(Exception):
    """Chemistry is not valid"""

    ...


class UnsupportedChemistryError(Exception):
    """Chemistry is not supported"""

    ...
