import csv

import go_obo_parser

import icdb


# Path for DB
PATH_DB = 'ic.db'

# Constants to be stored in the DB
ASSAY_IMMUNO = 'immuno'
ASSAY_MICROARRAY = 'microarray'
ASSAY_RNASEQ = 'RNA-seq'
DATASET_BIOGPS = 'U133AGNF1B'
DATASET_HPA = 'hpa'
EXTDB_BIOGPS = 'biogps'
EXTDB_BRENDA = 'brenda'
EXTDB_CHANNELPEDIA = 'channelpedia'
EXTDB_CHEMBL = 'chembl'
EXTDB_DRUGBANK = 'drugbank'
EXTDB_HMDB = 'hmdb'
EXTDB_HPA = 'hpa'
EXTDB_ZINC15 = 'zinc15'

# Paths for data files
PATH_BIOGPS_GCRMA = 'data/biogps/U133AGNF1B.gcrma.avg.csv'
PATH_BIOGPS_GNF1H_ANNOT = 'data/biogps/gnf1h.annot2007.tsv'
PATH_BIOGPS_TRANSLATION_TABLE = 'data/biogps/biogps_translation_table.csv'
PATH_BTO_BRENDA_TISSUE_OBO = 'data/bto/BrendaTissueOBO.obo'
PATH_BTO_TISSUES = 'data/bto/tissues.txt'
PATH_CHANNEL_CLASSES = 'data/channel-classes/channel-classes.csv'
PATH_CHANNELPEDIA_INFO = 'data/channel-classes/channelpedia_info.csv'
PATH_CHEMBL_ASSAYS_COMPOUND = 'data/chembl-assays/output_compound-assays_v1.1.dat'
PATH_CHEMBL_ASSAYS_DRUG = 'data/chembl-assays/output_drug-assays_v1.1.dat'
PATH_CHEMBL_HUMAN_ASSAYS_COMPOUND = 'data/chembl-assays/output_human_compound-assays_v1.1.dat' # Not used
PATH_CHEMBL_HUMAN_ASSAYS_DRUG = 'data/chembl-assays/output_human_drug-assays_v1.1.dat' # Not used
PATH_GO_ION_CHANNELS = 'data/go/ion-channels-1.csv'
PATH_GO_QUICKGO = 'data/go/QuickGO-ion-channel-COMBINED-human.dat'
PATH_HPA_NORMAL = 'data/hpa/normal_tissue.tsv'
PATH_HPA_PATHOLOGY = 'data/hpa/pathology.tsv'
PATH_HPA_RNA_TISSUE = 'data/hpa/rna_tissue.tsv'
PATH_HPA_TRANSLATION_TABLE = 'data/hpa/hpa_tissue_translation_table.csv'
PATH_TARGET_COMPOUND = 'data/target-compound/target-compound.csv' # Not used
PATH_TC_CHEMBL_COMPOUND = 'data/target-compound/output-query_ChEMBL-uniprot-compound.dat' # Not used
PATH_TC_CHEMBL_DRUG = 'data/target-compound/output-query_ChEMBL-uniprot-drug.dat' # Not used
PATH_TC_DRUG_INFO_DRUGBANK = 'data/target-compound/output-organized_drug-info-Drugbank.dat' # Not used
PATH_TC_DRUG_INFO_HMDB = 'data/target-compound/output-organized_drug-info-HMDB.dat' # Not used
PATH_TC_TARGET_INFO_DRUGBANK = 'data/target-compound/output-organized_target-info-Drugbank.dat' # Not used
PATH_TC_TARGET_INFO_HMDB = 'data/target-compound/output-organized_target-info-HMDB.dat' # Not used
PATH_TC_TTD = 'data/target-compound/output-organized-TTD.dat' # Not used
PATH_UNIPROT_ENSEMBL_IDMAPPING = 'data/uniprot/HUMAN_9606_idmapping_Ensembl.dat' # Not used
PATH_UNIPROT_GENBANK = 'data/uniprot/GenBankUniProt.txt' 
PATH_UNIPROT_HUMAN_IDMAPPING = 'data/uniprot/HUMAN_9606_idmapping.dat'


# DB setup routines

def setup_external_db_table():
    db = icdb.IonChannelDatabase(PATH_DB)
    db.create_external_db_table()
    db.add_external_db(EXTDB_BIOGPS, 'http://biogps.org')
    db.add_external_db(EXTDB_CHEMBL, 'https://www.ebi.ac.uk/chembl/')
    db.add_external_db(EXTDB_BRENDA, 'http://www.brenda-enzymes.org')
    db.add_external_db(EXTDB_DRUGBANK, 'https://www.drugbank.ca')
    db.add_external_db(EXTDB_ZINC15, 'http://zinc15.docking.org')
    db.add_external_db(EXTDB_HMDB, 'http://www.hmdb.ca')
    db.add_external_db(EXTDB_CHANNELPEDIA, 'http://channelpedia.epfl.ch')


def setup_tissue_table():
    db = icdb.IonChannelDatabase(PATH_DB)
    db.create_tissue_table()
    tissue_dict = {}
    tissue_file = open(PATH_BTO_TISSUES)
    for line in tissue_file:
        tissue_dict[line.strip()] = 1
    tissue_file.close()
    parser = go_obo_parser.parseGOOBO(PATH_BTO_BRENDA_TISSUE_OBO)
    for record in parser:
        tissue_name = record['name']
        if tissue_name in tissue_dict:
            bto_id = record['id']
            db.add_tissue(tissue_name, bto_id, '')


def input_biogps_tissue_data():
    db = icdb.IonChannelDatabase(PATH_DB)
    translation_file = open(PATH_BIOGPS_TRANSLATION_TABLE, 'rU')
    translation_reader = csv.reader(translation_file)
    translation_reader.next()
    for row in translation_reader:
        biogps_tissue = row[0]
        bto_tissue = row[1]
        db.add_db_tissue(EXTDB_BIOGPS, bto_tissue, biogps_tissue)
    translation_file.close()


def input_hpa_tissue_data():
    db = icdb.IonChannelDatabase(PATH_DB)
    translation_file = open(PATH_HPA_TRANSLATION_TABLE, 'rU')
    translation_reader = csv.reader(translation_file)
    translation_reader.next()
    for row in translation_reader:
        hpa_tissue = row[0]
        bto_tissue = row[1]
        db.add_db_tissue(EXTDB_HPA, bto_tissue, hpa_tissue)
    translation_file.close()


def setup_db_tissue_table():
    db = icdb.IonChannelDatabase(PATH_DB)
    db.create_db_tissue_table()
    input_biogps_tissue_data()
    input_hpa_tissue_data()


def setup_genbank_uniprot_table():
    db = icdb.IonChannelDatabase(PATH_DB)
    db.create_genbank_uniprot_table()
    gu_file = open(PATH_UNIPROT_GENBANK, 'rU')
    line_num = 0
    for line in gu_file:
        if line_num != 0:
            line_split = line.split('\t')
            if len(line_split) == 2:
                gbac = line_split[0].strip()
                upac = line_split[1].strip()
                if upac != '':
                    if not db.exists_genbank_uniprot(gbac, upac):
                        db.add_genbank_uniprot(gbac, upac)
        line_num += 1
    gu_file.close()


def input_chembl_compound_data():
    db = icdb.IonChannelDatabase(PATH_DB)
    compound_file = open(PATH_CHEMBL_ASSAYS_COMPOUND)
    for line in compound_file:
        line_split = line.split('\t')
        upac = line_split[0].strip()
        target_chembl_id = line_split[1].strip()
        target_type = line_split[2].strip()
        target_name = line_split[3].strip()
        compound_chembl_id = line_split[4].strip()
        compound_name = line_split[5].strip()
        iupac_name = line_split[6].strip()
        assay_chembl_id = line_split[7].strip()
        assay_standard_type = line_split[8].strip()
        assay_standard_relation = line_split[9].strip()
        assay_value = line_split[10].strip()
        assay_units = line_split[11].strip()
        assay_type = line_split[12].strip()
        assay_description = line_split[13].strip()
        if not db.exists_compound_by_chembl_id(compound_chembl_id):
            db.add_compound(
                compound_name, '', '', chembl_id=compound_chembl_id,
                synonyms=iupac_name, sourcedb_name=EXTDB_CHEMBL)
        if not db.exists_protein(upac):
            db.add_protein(
                upac, name=target_name, process_function=target_type,
                chembl_id=target_chembl_id)
        compound_id = db.get_compound_id_by_chembl_id(compound_chembl_id)
        db.add_interaction(
            upac, compound_id, action_type='',
            strength=assay_value.decode('utf_8').encode('utf_8'),
            strength_units=assay_units.decode('utf_8').encode('utf_8'),
            assay_type=assay_standard_type.decode('utf_8').encode('utf_8'),
            chembl_id=assay_chembl_id, sourcedb_name=EXTDB_CHEMBL)
    compound_file.close()


def input_chembl_drug_data():
    db = icdb.IonChannelDatabase(PATH_DB)
    drug_file = open(PATH_CHEMBL_ASSAYS_DRUG)
    for line in drug_file:
        line_split = line.split('\t')
        upac = line_split[0].strip()
        target_chembl_id = line_split[1].strip()
        target_type = line_split[2].strip()
        target_name = line_split[3].strip()
        drug_mechanism = line_split[4].strip()
        dosed_compound_chembl_id = line_split[5].strip()
        dosed_compound_name = line_split[6].strip()
        active_compound_chembl_id = line_split[7].strip()
        active_compound_name = line_split[8].strip()
        assay_chembl_id = line_split[9].strip()
        assay_standard_type = line_split[10].strip()
        assay_standard_relation = line_split[11].strip()
        assay_value = line_split[12].strip()
        assay_units = line_split[13].strip()
        assay_type = line_split[14].strip()
        assay_description = line_split[15].strip()
        if not db.exists_compound_by_chembl_id(active_compound_chembl_id):
            db.add_compound(
                '', '', active_compound_name,
                chembl_id=active_compound_chembl_id,
                synonyms=dosed_compound_name, sourcedb_name=EXTDB_CHEMBL)
        if not db.exists_protein(upac):
            db.add_protein(
                upac, name=target_name, process_function=target_type,
                chembl_id=target_chembl_id)
        compound_id = db.get_compound_id_by_chembl_id(
            active_compound_chembl_id)
        db.add_interaction(
            upac, compound_id,
            action_type=drug_mechanism.decode('utf_8').encode('utf_8'),
            strength=assay_value.decode('utf_8').encode('utf_8'),
            strength_units=assay_units.decode('utf_8').encode('utf_8'),
            assay_type=assay_standard_type.decode('utf_8').encode('utf_8'),
            chembl_id=assay_chembl_id, sourcedb_name=EXTDB_CHEMBL)
    drug_file.close()


def input_chembl_assay_data():
    input_chembl_compound_data()
    input_chembl_drug_data()


def input_ion_channel_data():
    db = icdb.IonChannelDatabase(PATH_DB)
    prot_file = open(PATH_GO_ION_CHANNELS, 'rU')
    prot_reader = csv.reader(prot_file)
    prot_reader.next()
    for row in prot_reader:
        protein_name = row[2]
        process_function = row[3]
        upac = row[0]
        if upac != '':
            gene_symbol = row[1]
            if gene_symbol == '-':
                gene_symbol = ''
            ions = row[4]
            gating = row[5]
            ion_channel_subclass = row[6]
            in_betse = row[7]
            if in_betse in ['yes', 'soon', 'doable']:
                in_betse = 'Y'
            else:
                in_betse = 'N'
            if not db.exists_protein(upac):
                db.add_protein(
                    upac, gene_symbol, protein_name, process_function, ions,
                    gating, in_betse, ion_channel_subclass)
            else:
                old_protein_record = db.lookup_protein(upac)
                old_upac = old_protein_record[0]
                old_gene_symbol = old_protein_record[1]
                old_name = old_protein_record[2]
                old_process_function = old_protein_record[3]
                old_ions = old_protein_record[4]
                old_gating = old_protein_record[5]
                old_in_betse = old_protein_record[6]
                old_ion_channel_class_desc = old_protein_record[7]
                old_ion_channel_subclass = old_protein_record[8]
                old_chembl_id = old_protein_record[9]
                db.update_protein_gene_symbol(upac, gene_symbol)
                db.update_protein_name(upac, protein_name)
                db.update_protein_ions(upac, ions)
                db.update_protein_gating(upac, gating)
                if old_in_betse == '':
                    db.update_protein_in_betse(upac, in_betse)
                db.update_protein_ion_channel_sub_class(
                    upac, ion_channel_subclass)
    prot_file.close()


def setup_protein_compound_interaction_tables():
    db = icdb.IonChannelDatabase(PATH_DB)
    db.create_protein_table()
    db.create_compound_table()
    db.create_interaction_table()
    input_chembl_assay_data()
    input_ion_channel_data()


def get_biogps_probeset_id_to_gene_symbol_dict():
    annot_dict = {}
    chip_annot_file = open(PATH_BIOGPS_GNF1H_ANNOT, 'rU')
    chip_annot_reader = csv.reader(chip_annot_file, delimiter='\t')
    chip_annot_reader.next()
    for row in chip_annot_reader:
        probeset_id = row[0]
        gene_symbol = row[6]
        if gene_symbol != '':
            annot_dict[probeset_id] = gene_symbol
    chip_annot_file.close()
    return annot_dict


def input_biogps_expression_data():
    db = icdb.IonChannelDatabase(PATH_DB)
    annot_dict = get_biogps_probeset_id_to_gene_symbol_dict()
    microarray_file = open(PATH_BIOGPS_GCRMA, 'rU')
    microarray_reader = csv.reader(microarray_file)
    header = microarray_reader.next()
    tissue_list = header[1:]
    for row in microarray_reader:
        probeset_id = row[0]
        if probeset_id in annot_dict:
            gene_symbol = annot_dict[probeset_id]
            upac_resultset = db.get_uniprot_accnums_by_gene_symbol(
                gene_symbol.strip())
            if len(upac_resultset) != 0:
                for i in xrange(0, len(tissue_list)):
                    expr_level = float(row[i + 1])
                    tissue_name = db.get_tissue_name_by_db_tissue_name(
                        EXTDB_BIOGPS, tissue_list[i])
                    for upac_result in upac_resultset:
                        upac = upac_result[0]
                        db.add_expression(
                            tissue_name, upac, expr_level,
                            assay_type=ASSAY_MICROARRAY,
                            dataset_name=DATASET_BIOGPS, 
                            sourcedb_name=EXTDB_BIOGPS)
    microarray_file.close()


# Function is not used
def get_ensembl_to_uniprot_dict():
    ensembl_to_uniprot = {}
    map_file = open(PATH_UNIPROT_ENSEMBL_IDMAPPING)
    for line in map_file:
        row = line.strip().split('\t')
        upac = row[0]
        ensg = row[2]
        ensembl_to_uniprot[ensg] = upac
    map_file.close()
    return ensembl_to_uniprot


def input_hpa_rna_expression_data():
    db = icdb.IonChannelDatabase(PATH_DB)
    hpa_file = open(PATH_HPA_RNA_TISSUE)
    hpa_file.next()
    for line in hpa_file:
        row = line.strip().split('\t')
        if len(row) >= 5 and row[0].strip() != '':
            gene_symbol = row[1]
            upac_record_list = db.get_uniprot_accnums_by_gene_symbol(
                gene_symbol)
            if len(upac_record_list) != 0:
                for upac_record in upac_record_list:
                    upac = upac_record[0]
                    gene_name = row[1]
                    tissue_name = row[2]
                    tissue_bto = db.get_tissue_name_by_db_tissue_name(
                        EXTDB_HPA, tissue_name)
                    expr_value = float(row[3])
                    db.add_expression(
                        tissue_bto, upac, expr_value,
                        assay_type=ASSAY_RNASEQ, dataset_name=DATASET_HPA,
                        sourcedb_name=EXTDB_HPA)
    hpa_file.close()


def input_hpa_cancer_expression_data():
    db = icdb.IonChannelDatabase(PATH_DB)
    hpa_file = open(PATH_HPA_PATHOLOGY)
    hpa_file.next()
    expr_qual_choices = ['High', 'Medium', 'Low', 'Not detected']
    for line in hpa_file:
        row = line.strip().split('\t')
        if len(row) >= 6:
            gene_symbol = row[1]
            upac_record_list = db.get_uniprot_accnums_by_gene_symbol(
                gene_symbol)
            if len(upac_record_list) != 0:
                for upac_record in upac_record_list:
                    upac = upac_record[0]
                    tissue_name = row[2]
                    tissue_bto = db.get_tissue_name_by_db_tissue_name(
                        EXTDB_HPA, tissue_name)
                    expr_value = 0.0
                    data_is_okay = False
                    try:
                         patient_counts = map(int, row[3:7])
                         data_is_okay = True
                    except:
                        pass
                    if data_is_okay:
                        patient_count_max = max(patient_counts)
                        patient_count_max_idx = [
                            idx for idx, value in enumerate(patient_counts)
                            if value==patient_count_max][-1]
                        expr_value_qual = expr_qual_choices[
                            patient_count_max_idx]
                        db.add_expression(
                            tissue_bto, upac, expr_value,
                            expr_level_qual=expr_value_qual,
                            assay_type=ASSAY_IMMUNO, dataset_name=DATASET_HPA,
                            sourcedb_name=EXTDB_HPA)
    hpa_file.close()


def input_hpa_normal_expression_data():
    db = icdb.IonChannelDatabase(PATH_DB)
    hpa_file = open(PATH_HPA_NORMAL)
    hpa_file.next()
    for line in hpa_file:
        row = line.strip().split('\t')
        if len(row) >= 6:
            gene_symbol = row[1]
            upac_record_list = db.get_uniprot_accnums_by_gene_symbol(
                gene_symbol)
            if len(upac_record_list) == 6:
                for upac_record in upac_record_list:
                    upac = upac_record[0]
                    tissue_name = row[2]
                    tissue_bto = db.get_tissue_name_by_db_tissue_name(
                        EXTDB_HPA, tissue_name)
                    expr_value = 0.0
                    expr_value_qual = row[4]
                    if expr_value_qual in [
                            'Not detected', 'Low', 'Medium', 'High']:
                        db.add_expression(
                            tissue_bto, upac, expr_value,
                            expr_level_qual=expr_value_qual,
                            assay_type=ASSAY_IMMUNO, dataset_name=DATASET_HPA,
                            sourcedb_name=EXTDB_HPA)
    hpa_file.close()


def setup_expression_table():
    db = icdb.IonChannelDatabase(PATH_DB)
    db.create_expression_table()
    input_biogps_expression_data()
    input_hpa_rna_expression_data()
    input_hpa_cancer_expression_data()
    input_hpa_normal_expression_data()


def setup_channel_class_tables():
    db = icdb.IonChannelDatabase(PATH_DB)
    db.create_channel_super_class_table()
    db.create_channel_class_table()
    db.create_channel_sub_class_table()
    channelpedia_dict = {}
    channelpedia_file = open(PATH_CHANNELPEDIA_INFO)
    channelpedia_reader = csv.reader(channelpedia_file)
    channelpedia_reader.next()
    for row in channelpedia_reader:
        subclass = row[1]
        intro_text = row[2]
        url = row[3]
        channelpedia_dict[subclass] = (intro_text, url)
    channelpedia_file.close()
    class_file = open(PATH_CHANNEL_CLASSES, 'rU')
    class_reader = csv.reader(class_file)
    class_reader.next()
    for row in class_reader:
        channel_superclass = row[0]
        channel_class = row[1]
        channel_subclass = row[2]
        subfamily = row[3]
        gene_symbol = row[4]
        if not db.exists_channel_super_class(channel_superclass):
            db.add_channel_super_class(channel_superclass)
        if not db.exists_channel_class(channel_class):
            db.add_channel_class(channel_class, channel_superclass)
        if not db.exists_channel_sub_class(channel_subclass):
            if channel_subclass in channelpedia_dict:
                channelpedia_info, channelpedia_url = channelpedia_dict[
                    channel_subclass]
            else:
                channelpedia_info = ''
                channelpedia_url = ''
            db.add_channel_sub_class(
                channel_subclass, channel_class,
                channelpedia_text=channelpedia_info.decode('utf_8').encode('utf_8'),
                channelpedia_url=channelpedia_url.decode('utf_8').encode('utf_8'))
        upac_record_list = db.get_uniprot_accnums_by_gene_symbol(gene_symbol)
        for upac_record in upac_record_list:
            upac = upac_record[0]
            db.update_protein_sub_class(upac, channel_subclass)
    class_file.close()


def setup_go_term_table():
    db = icdb.IonChannelDatabase(PATH_DB)
    db.create_go_term_table()
    go_file = open(PATH_GO_QUICKGO)
    go_reader = csv.reader(go_file, delimiter='\t')
    go_reader.next()
    for row in go_reader:
        db_name = row[0]
        if db_name == 'UniProtKB':
            upac = row[1]
            goid = row[4]
            if db.exists_protein(upac):
                db.add_go_term(upac, goid)
    go_file.close()


def setup_specificity_table_with_jp_method():
    db = icdb.IonChannelDatabase(PATH_DB)
    db.create_specificity_table()
    all_upac_list = []
    upac_record_list = db.get_protein_uniprot_accnums()
    for upac_record in upac_record_list:
        upac = upac_record[0]
        all_upac_list.append(upac)
    tissue_name_record_list = db.get_tissue_names()
    all_tissue_name_list = []
    for tissue_name_record in tissue_name_record_list:
        tissue_name = tissue_name_record[0]
        all_tissue_name_list.append(tissue_name)
    specificity_score_threshold = 5.0
    specificity_table = {}  
    expr_table = {}
    for upac in all_upac_list:
        specificity_score = 0
        ntissues = 0
        for tissue_name in all_tissue_name_list:
            expr_level_record_list = \
                db.get_expression_level_by_uniprot_accnum_tissue_dataset(
                upac, tissue_name, EXTDB_HPA)
            if expr_level_record_list:
                expr_level_rna = expr_level_record_list[0][0]
                expr_table[(tissue_name, upac)] = expr_level_rna
                if expr_level_rna > specificity_score_threshold:
                    specificity_score +=1
                ntissues += 1
        if ntissues > 0:
            specificity_score = 1 - specificity_score*1./ntissues
        else:
            specificity_score = None
        specificity_table[upac] = specificity_score
    for tissue_name, upac in expr_table:
        expr_level_rna = expr_table[(tissue_name, upac)]
        specificity_score_tissue = specificity_table[upac] * expr_level_rna
        if specificity_score_tissue:
            db.add_specificity(tissue_name, upac, specificity_score_tissue)


def vacuum_db():
    db = icdb.IonChannelDatabase(PATH_DB)
    db.vacuum_db()


def dump_db():
    db = icdb.IonChannelDatabase(PATH_DB)
    db.dump_db()


def print_db_stats():
    db = icdb.IonChannelDatabase(PATH_DB)
    db.print_db_stats()
