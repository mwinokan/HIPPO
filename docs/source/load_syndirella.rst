
==========================
Loading Syndirella outputs
==========================

Syndirella has been developed to produce a machine/program friendly output:

1. `final_output.csv` which summarises all the elaboration series per scaffold

2. HIPPO-friendly pickled dataframes for each scaffold in the syntax:

::

    {inchikey}_{reaction_uuid}_to_hippo.pkl.gz

Parsing the summary csv
=======================

To get the reference structures construct a mapping between the pickle paths and the shortcode:

::

    import pandas as pd
    import numpy as np

    df = pd.read_csv("/path/to/final_output.csv")

    ref_lookup = {}

    for i, row in df.iterrows():
        to_hippo, template = row["to_hippo"], row["template"]
        if to_hippo is np.nan:
            continue

        # reconstruct the shortcode from the longcode
        ref_lookup[to_hippo] = Path(template).name.split("_")[0].replace("x", "")

Loading the base/scaffold compounds from each pickle
====================================================

::

    for to_hippo, template in ref_lookup.items():
        df = animal.add_syndirella_elabs(to_hippo, base_only=True, check_chemistry=True, reference=animal.poses[template])


Loading all the elaboration series
==================================

Due to the large number of molecules this is best submitted as a batch job on a cluster. Example python script:

::

    import hippo
    import pandas as pd
    from pathlib import Path
    import os
    import numpy as np

    from mlog import setup_logger
    logger = setup_logger('CHIKV_prod2')

    os.system("cp -v CHIKV_prod1.sqlite CHIKV_prod2.sqlite")

    animal = hippo.HIPPO("CHIKV_prod2", "CHIKV_prod2.sqlite")

    df = pd.read_csv(
        "/opt/xchem-fragalysis-2/kfieseler/CHIKV-Mac-syndirella-run/sept9_syndirella_final_output.csv"
    )

    ref_lookup = {}

    for i, row in df.iterrows():
        to_hippo, template = row["to_hippo"], row["template"]
        if to_hippo is np.nan:
            continue
        ref_lookup[to_hippo] = Path(template).name.split("_")[0].replace("x", "")

    for i, (to_hippo, template) in enumerate(ref_lookup.items()):
        logger.title(i)
        logger.reading(to_hippo)
        df = animal.add_syndirella_elabs(to_hippo, reference=animal.poses[template])

    animal.db.close()
