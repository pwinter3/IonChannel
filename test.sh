#!/bin/bash

set -x

rm levin.db
python -c "from levindb import *;create_tables()"
python -c "from levindb import *;cleanup_externaldb()"
python -c "from levindb import *;setup_externaldb()"
python -c "from levindb import *;cleanup_ontology()"
python -c "from levindb import *;setup_ontology()"
python -c "from levindb import *;cleanup_dbtissue()"
python -c "from levindb import *;setup_dbtissue()"
python -c "from levindb import *;cleanup_genbank_uniprot()"
python -c "from levindb import *;setup_genbank_uniprot()"
python -c "from levindb import *;cleanup_proteins()"
python -c "from levindb import *;cleanup_compounds()"
python -c "from levindb import *;cleanup_interactions()"
python -c "from levindb import *;input_chembl_assays()"
python -c "from levindb import *;setup_ion_channels()"
python -c "from levindb import *;cleanup_expression()"
python -c "from levindb import *;process_biogps_data()"
python -c "from levindb import *;process_hpa_rna_data()"
python -c "from levindb import *;process_hpa_cancer_data()"
python -c "from levindb import *;process_hpa_normal_data()"
python -c "from levindb import *;cleanup_channel_classes()"
python -c "from levindb import *;setup_channel_classes()"
python -c "from levindb import *;cleanup_specificity()"
python -c "from levindb import *;calculate_specificity_jp()"
python -c "from levindb import *;dump_database()"
python -c "from levindb import *;print_db_stats()"
python -c "from levindb import *;vacuum()"

set +x
