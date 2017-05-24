#!/bin/bash
#rm levin.db
python -c "from levindb import *;create_tables()"
python -c "from levindb import *;cleanup_externaldb()"
python -c "from levindb import *;setup_externaldb()"
python -c "from levindb import *;cleanup_ontology()"
python -c "from levindb import *;setup_ontology()"
python -c "from levindb import *;cleanup_proteins()"
python -c "from levindb import *;setup_proteins()"
python -c "from levindb import *;cleanup_dbtissue()"
python -c "from levindb import *;setup_biogps()"
python -c "from levindb import *;cleanup_expression()"
python -c "from levindb import *;process_biogps_data()"
python -c "from levindb import *;cleanup_compounds()"
python -c "from levindb import *;setup_chembl()"
python -c "from levindb import *;dump_database()"
