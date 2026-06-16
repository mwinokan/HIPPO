## Testing new HIPPO

Hippo testing during active development phase: mount the code in git branch directly to container and have it handle notebook server.


## Environment variables
Not strictly necessary but it's convenient to have database connection parameters in `.env` file, in project root directory. Should contain minimally:

```
DB_NAME=designdb
DB_USER=postgres
DB_PASSWORD=s_URzt7CWfWZ.AXD7RcF
DB_HOST=database

TA_AUTH_SERVICE=https://ta-authenticator.xchem.diamond.ac.uk/
TA_AUTH_QUERY_KEY=<auth key>
```

These are for local connection, values may be different for kubernetes deployment


## Running services

NB! depending on your environment, you may have to prefix your docker commands with `sudo`

### Build the app container
In project root directory:

```bash
docker build --no-cache . -t hippo_backend:latest
```

`--no-cache` may not be necessary, but I find that without that it sometimes does not upgrade numpy.


### If using local database, build the designdb postgres container
In `<project root>/images/xchem-designdb` directory:

```bash
docker build -t xchem_designdb:latest .
```
This may take some time


### Launch the service(s)

`docker-compose.yaml` contains instructions for docker to run the services.

```bash
docker compose up
```

If not using local database, just the app container:

```bash
docker compose up backend
```

After running the `up` command, at the very end you should see notebook server addres:

```

hippo_backend   |     Or copy and paste one of these URLs:
hippo_backend   |         http://localhost:8888/lab?token=90de0d2f297079ee393cdc202c06064406c5cb3a8032e8d8
hippo_backend   |         http://127.0.0.1:8888/lab?token=90de0d2f297079ee393cdc202c06064406c5cb3a8032e8d8
hippo_backend   | [I 2026-04-23 09:49:32.420 ServerApp] Skipped non-installed server(s): basedpyright, bash-language-server, dockerfile-language-server-nodejs, javascript-typescript-langserver, jedi-language-server, julia-language-server, pyrefly, pyright, python-language-server, python-lsp-server, r-languageserver, sql-language-server, texlab, typescript-language-server, unified-language-server, vscode-css-languageserver-bin, vscode-html-languageserver-bin, vscode-json-languageserver-bin, yaml-language-server

```
Copy one of the addresses to a broser tab and you're in jupyter lab environment.


### Cleaning up when done

```bash
docker compose down
```

When using local database and it's necessary to wipe db contents:

```bash
docker compose down -v
```


## Test commands you originally shared with me

I tried to keep the commands as they were before but there are some changes, and considering the upcoming permissions ands scope issues, there will be more.

### Imports

Original imports

```
import hippo
```

New imports

```
from hippo import HIPPO
```
(+ whatever else you might need)



### Setup animal
Original setup

```
target_name = "Flavi_NS5_RdRp"
animal = hippo.HIPPO(target_name, f"{target_name}.sqlite", update_legacy=True)
```

New animal setup:

```
target_name = "Flavi_NS5_RdRp"
target_access_string = "lb18145-1"
username = '<your FedID>'

animal = hippo.HIPPO(
    target_name=target_name,
    target_access_string=target_access_string,
    username=username,
)


```

If you didn't set up an `.env` file, you need to give connection parameters here (using the appropriate values of course):

```
db = {
    'DB_NAME': designdb,
    'DB_USER': postgres,
    'DB_PASSWORD': s_URzt7CWfWZ.AXD7RcF,
    'DB_HOST': database,
    'POSTGRES_PORT': '5432',
}
animal = hippo.HIPPO(
    target_name=target_name,
    target_access_string=target_access_string,
    username=username,
    db=db,
)
```

Unlike before, `HIPPO` is now a bootstrap function that sets up database. Since we're using it as a library not a standalone app, it needs to be called by user. That means, many if not most library objects cannot be imported before the animal is initialised. For example `RouteSet` (used below), can be imported only now:

```
from designdb.sets.route import RouteSet
```



Most important change here (even though it's not visible) - the `target_name` argument, which used to be just a project name, now actually is a target name. If there's no target by that name, it will be created, otherwise it will be fetched from the db. `animal` will remain associated with this target during it's lifecycle, you cannot work with other targets once initialised.


`update_legacy` flag doesn't work at the moment and going forward, probably won't be necessary at all.

It's still possible to work with local sqlite databases. Let me know if you're interested in this, I may have to do some tweaks to make it more convenient.


### Registering methods

The following methods seem necessary for the test workflow
```
animal.register_enumeration_method(name="fragmenstein", version="1.0.0", description="XChem's implementation of Fragmenstein merges")


animal.register_pose_method(name="fragmenstein", version="1.0.0", description="XChem's implementation of Fragmenstein placement")
animal.register_pose_method(name="xray", version="1.0.0", description="Appropriate description")
animal.register_pose_method(name="gnina_repose", version="1.0.0", description="Appropriate description")


animal.register_scoring_method(name="gnina_cnn_vs", version="1.3.2", description="GNINA CNN VS score")
animal.register_scoring_method(name="fragmenstein_energy", version="1.0.0", description="appropriate description")
animal.register_scoring_method(name="fragmenstein_distance", version="1.0.0", description="appropriate description")
animal.register_scoring_method(name="moc_combo_multiref", version="1.3.2", description="Symmetric similarity score (equal weighting given to references and queries) with given reference/inspiration compounds set as references, and placed designs as the query")
animal.register_scoring_method(name="moc_combo_multiref", version="0.1.0", description="Symmetric similarity score (equal weighting given to references and queries) with given reference/inspiration compounds set as references, and placed designs as the query")
```


### Add fragalysis hits
Original command:

```
animal.add_hits(
    target_name=target_name,
    metadata_csv=f"{target_name}/metadata.csv",
    aligned_directory= f"{target_name}/aligned_files",
    load_pose_mols=True,
)
```

New command:

```
animal.add_hits(
    metadata_csv=f"{target_name}/metadata.csv",
    aligned_directory= f"{target_name}/aligned_files",
)
```

Since the animal already knows about the target, there' no need to specify it.
`load_pose_mols` is gone because pose registration is done quite differently now.


## Create input for Fragmenstein and Knitwork
No changes here:

```
fragment_hits = animal.poses(tag="hits")
fragment_hits

fragment_hits.write_sdf("fragment_hits.sdf")
fragment_hits.to_knitwork("knitwork_input.csv", aligned_files_dir="aligned_files")
```

### Load BulkDock SDF
Original code:

```
SDFs = [
    "openbind_flavi_ns5_rdrp_c1_fragmenstein_split2000_batch002_820954.sdf",
    "openbind_flavi_ns5_rdrp_c1_fragmenstein_split2000_batch000_820952.sdf",
    "openbind_flavi_ns5_rdrp_c1_fragmenstein_split2000_batch001_820953.sdf"
]

for sd in SDFs:
    full_path = os.path.join("bulkdock", sd)

    key = sd.strip(".sdf")

    animal.load_sdf(
        target= target_name,
        path=full_path,
        inspiration_col="inspiration_ids",
        reference_col="reference_id",
        compound_tags = [key],
        pose_tags = ["fragmenstein_placed", key],
        name_col = "ID",
    )
```

New code, `target_name` is gone:

```
SDFs = [
    "openbind_flavi_ns5_rdrp_c1_fragmenstein_split2000_batch002_820954.sdf",
    "openbind_flavi_ns5_rdrp_c1_fragmenstein_split2000_batch000_820952.sdf",
    "openbind_flavi_ns5_rdrp_c1_fragmenstein_split2000_batch001_820953.sdf"
]

for sd in SDFs:
    full_path = os.path.join("bulkdock", sd)

    key = sd.strip(".sdf")

    animal.load_sdf(
        path=full_path,
        inspiration_col="inspiration_ids",
        reference_col="reference_id",
        compound_tags = [key],
        pose_tags = ["fragmenstein_placed", key],
        name_col = "ID",
    )
```

### Generate Fragalysis RHS input
No changes

```
poses = animal.poses.get_by_tag("fragmenstein_placed")

poses.to_fragalysis("flavi_ns5_rdrp_bulkdock_poses.sdf", method="fragmenstein", submitter_name = "Lauren Reid", submitter_email= "lauren.reid@medchemica.com", submitter_institution="MedChemica", copy_reference_pdbs=True)
```


### Load GNINA poses and scores
Drops `target_name`, otherwise no changes

```
animal.load_sdf(
    path="gnina/output_sdfs/z0625b_ligands_minimized.sdf",
    pose_tags = ["gnina_minimised"],
)
```


### Create Syndirella inputs
No changes

```
poses.to_syndirella("syndirella_input.csv")
```

### Load Syndirella retrosynthesis routes
No changes

```
animal.add_syndirella_routes(
    "syndirella/retro/justretroquery_manifold_ZWWZBAUAVBEARZ-UHFFFAOYSA-N-scaffold-A.pkl.gz",
    CAR_only=False,
    check_chemistry=False,
)
```


### Load Syndirella elaborations
Old code:

```
comp = animal.compounds.get_by_smiles("CN(C[C@H]1CCNC1)c1ccc2ccccc2n1")

print(comp.id)
print(comp.smiles)

routes = hippo.RouteSet.from_product_ids(animal.db, [comp.id])

assert len(routes) == 1, "Wrong number of routes"

route = routes.pop()

animal.add_syndirella_elabs("syndirella/elabs/ZWWZBAUAVBEARZ-UHFFFAOYSA-N_a7d7696daae73aca44078d04fc8c3093_structured_output.pkl.gz", scaffold_route=route)
```

New code:

```
comp = animal.compounds.get_by_smiles("CN(C[C@H]1CCNC1)c1ccc2ccccc2n1")

print(comp.id)
print(comp.smiles)

routes = RouteSet.from_product_ids(animal.db, [comp.id])

assert len(routes) == 1, "Wrong number of routes"

route = routes.pop()

animal.add_syndirella_elabs("syndirella/elabs/ZWWZBAUAVBEARZ-UHFFFAOYSA-N_a7d7696daae73aca44078d04fc8c3093_structured_output.pkl.gz", scaffold_route=route)
```

RouteSet must be imported and called directly, not through animal. This is something that probably needs to change, but I'm not sure which way. I'll know better once work on target and TAS scope begins.

NB! the last command did not run successfully for me, I didn't have original SDF files on disk, so it couldn't find the reference files, and even if I disabled that, there was nan reference in the data frame. If the files were indeed correct, I'll have to revisit that.
