#!/bin/bash

set -x

#rm levin.db
#python -c "from levindb import *;create_tables()"
#python -c "from levindb import *;cleanup_externaldb()"
#python -c "from levindb import *;setup_externaldb()"
#python -c "from levindb import *;cleanup_ontology()"
#python -c "from levindb import *;setup_ontology()"
#python -c "from levindb import *;cleanup_proteins()"
#python -c "from levindb import *;cleanup_dbtissue()"
#python -c "from levindb import *;setup_biogps()"
#python -c "from levindb import *;cleanup_compounds()"
#python -c "from levindb import *;cleanup_interactions()"
#python -c "from levindb import *;input_chembl_assays()"
#python -c "from levindb import *;setup_ion_channels()"
#python -c "from levindb import *;cleanup_expression()"
#python -c "from levindb import *;process_biogps_data()"
#python -c "from levindb import *;process_protein_atlas_data()"
#python -c "from levindb import *;setup_channel_classes()"
#python -c "from levindb import *;dump_database()"
python -c "from levindb import *;write_db_stats()"

# Not implemented or deprecated:
#python -c "from levindb import *;setup_target_compound()"
#python -c "from levindb import *;input_drugbank()"
#python -c "from levindb import *;input_hmdb()"
#python -c "from levindb import *;input_ttd()"
#python -c "from levindb import *;input_chembl_uniprot()"
#python -c "from levindb import *;setup_chembl()"
#python -c "from levindb import *;setup_all_proteins()"
#python -c "from levindb import *;get_expression_prob_for_all_tissues()"
#python -c "from levindb import *;write_important_data()"

set +x
