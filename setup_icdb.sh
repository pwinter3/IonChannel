#!/bin/bash
set -x
[ -e ic.db ] && rm ic.db
python -c "from setup_icdb import *;setup_external_db_table()"
python -c "from setup_icdb import *;setup_tissue_table()"
python -c "from setup_icdb import *;setup_db_tissue_table()"
python -c "from setup_icdb import *;setup_protein_compound_interaction_tables()" # Slow
python -c "from setup_icdb import *;input_in_betse_data()"
python -c "from setup_icdb import *;setup_go_term_table()"
python -c "from setup_icdb import *;setup_expression_table()" # Slow
python -c "from setup_icdb import *;setup_channel_class_tables()"
python -c "from setup_icdb import *;setup_specificity_table_with_jp_method()" # Slow
python -c "from setup_icdb import *;vacuum_db()"
python -c "from setup_icdb import *;dump_db()"
python -c "from setup_icdb import *;print_db_stats()"
set +x
