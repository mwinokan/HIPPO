
import pandas as pd
import molparse as mp
import mout

'''

Syndirella file structure on IRIS:

Home_dir: kfieseler/A71EV2A

|-- /routes_data: Contains master csvs of routes to compounds
|-- /elabs: Contains all the placements
    |-- /1_step_dec6_batched: 1 step cmpds started to place on DEC 6 (most number of base compounds)
        |-- /batch_1
            |-- /submitter_method_base_cmpd_id
                |-- submitter_method_base_cmpd_id_num_analogs_1_of_1.csv: master csv with all data of all elabs
                |-- submitter_method_base_cmpd_id_num_analogs_1_of_1_batch_[num_atom_low]-[num_atom_high].csv: data of elabs with specified num atom diff range
                |-- /success
                    |-- /cmpd_name
                    |-- success.csv: should be original csv with fragmenstein information parsed directly from json.
         % Continue rest of batches
    |-- /1_step_dec7_batched: 1 step cmpds started to place on DEC 6 (same dir structure)
    |-- /2_step: All 2 step cmpds (same did structure)

'''

def parse_routes_data(path):
	...

def parse_elabs_csv(path):

	df = pd.read_csv(path)

	

	compound = Compound.from_mol(comp_name, mol, tags=tags)

