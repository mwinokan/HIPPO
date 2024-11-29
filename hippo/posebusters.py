import mrich


mrich.debug("importing posebusters.PoseBusters")
from posebusters import PoseBusters

# setup posebusters
mrich.debug("Setting up PoseBusters(config='mol')")
MOL_ENGINE = PoseBusters(config="mol")
MOL_ENGINE._initialize_modules()
MOL_ENGINE.module_name.pop(0)
MOL_ENGINE.module_func.pop(0)
MOL_ENGINE.module_args.pop(0)
MOL_ENGINE.fname.pop(0)
MOL_ENGINE.module_test_names = {
    c["name"]: c["chosen_binary_test_output"] for c in MOL_ENGINE.config["modules"]
}


def check_mol(mol, debug: bool = False):

    failed = False

    for name, fname, func, args in zip(
        MOL_ENGINE.module_name,
        MOL_ENGINE.fname,
        MOL_ENGINE.module_func,
        MOL_ENGINE.module_args,
    ):

        assert len(args) == 1 and "mol_pred" in args

        results = func(mol_pred=mol)

        test_names = MOL_ENGINE.module_test_names[name]

        for test_name in test_names:

            passed = results["results"][test_name]

            if debug:
                mrich.var(f"{name}.{test_name}", results["results"][test_name])

            if not passed:
                failed = True
                mrich.error(f"Molecule failed {name}.{test_name}")

    return not failed
