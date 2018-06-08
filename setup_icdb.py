import csv

import go_obo_parser

import icdb


# Path for DB
PATH_DB = 'ic.db'

# Constants to be stored in the DB
ACTION_TYPE_ACTIVATOR = 'ACTIVATOR'
ACTION_TYPE_AGONIST = 'AGONIST'
ACTION_TYPE_ALLOSTERIC_ANTAGONIST = 'ALLOSTERIC ANTAGONIST'
ACTION_TYPE_ANTAGONIST = 'ANTAGONIST'
ACTION_TYPE_BLOCKER = 'BLOCKER'
ACTION_TYPE_INHIBITOR = 'INHIBITOR'
ACTION_TYPE_MODULATOR = 'MODULATOR'
ACTION_TYPE_NEGATIVE_ALLOSTERIC_MODULATOR = 'NEGATIVE ALLOSTERIC MODULATOR'
ACTION_TYPE_OPENER = 'OPENER'
ACTION_TYPE_PARTIAL_AGONIST = 'PARTIAL AGONIST'
ACTION_TYPE_POSITIVE_ALLOSTERIC_MODULATOR = 'POSITIVE ALLOSTERIC MODULATOR'
ACTION_TYPE_POSITIVE_MODULATOR = 'POSITIVE MODULATOR'
ACTION_TYPE_STABILISER = 'STABILISER'
EXPR_ASSAY_IMMUNO = 'immuno'
EXPR_ASSAY_MICROARRAY = 'microarray'
EXPR_ASSAY_RNASEQ = 'RNA-seq'
ASSAY_TYPE_ADMET = 'A'
ASSAY_TYPE_BINDING = 'B'
ASSAY_TYPE_FUNCTIONAL = 'F'
ASSAY_TYPE_TOXICITY = 'T'
ASSAY_TYPE_UNASSIGNED = 'U'
COMPOUND_TYPE_COMPOUND = 'compound'
COMPOUND_TYPE_DRUG = 'drug'
DATASET_BIOGPS = 'U133AGNF1B'
DATASET_HPA = 'hpa'
EXPR_H = 'High'
EXPR_L = 'Low'
EXPR_M = 'Medium'
EXPR_ND = 'Not detected'
EXTDB_BIOGPS = 'biogps'
EXTDB_BRENDA = 'brenda'
EXTDB_CHANNELPEDIA = 'channelpedia'
EXTDB_CHEMBL = 'chembl'
EXTDB_DRUGBANK = 'drugbank'
EXTDB_HMDB = 'hmdb'
EXTDB_HPA = 'hpa'
EXTDB_ZINC15 = 'zinc15'
NO = 'N'
YES = 'Y'
ACTION_TYPE_CHOICES = [
    ACTION_TYPE_ACTIVATOR,
    ACTION_TYPE_AGONIST,
    ACTION_TYPE_ALLOSTERIC_ANTAGONIST,
    ACTION_TYPE_ANTAGONIST,
    ACTION_TYPE_BLOCKER,
    ACTION_TYPE_INHIBITOR,
    ACTION_TYPE_MODULATOR,
    ACTION_TYPE_NEGATIVE_ALLOSTERIC_MODULATOR,
    ACTION_TYPE_OPENER,
    ACTION_TYPE_PARTIAL_AGONIST,
    ACTION_TYPE_POSITIVE_ALLOSTERIC_MODULATOR,
    ACTION_TYPE_POSITIVE_MODULATOR,
    ACTION_TYPE_STABILISER,
    ]
ASSAY_TYPE_CHOICES = [
    ASSAY_TYPE_ADMET,
    ASSAY_TYPE_BINDING,
    ASSAY_TYPE_FUNCTIONAL,
    ASSAY_TYPE_TOXICITY,
    ASSAY_TYPE_UNASSIGNED,
    ]

# Paths for data files
PATH_BIOGPS_GCRMA = 'data/biogps/U133AGNF1B.gcrma.avg.csv'
PATH_BIOGPS_GNF1H_ANNOT = 'data/biogps/gnf1h.annot2007.tsv'
PATH_BIOGPS_TRANSLATION_TABLE = 'data/biogps/biogps_translation_table.csv'
PATH_BTO_BRENDA_TISSUE_OBO = 'data/bto/BrendaTissueOBO.obo'
PATH_BTO_TISSUES = 'data/bto/tissues.txt'
PATH_CHANNEL_CLASSES = 'data/channel-classes/channel-classes.csv'
PATH_CHANNELPEDIA_INFO = 'data/channel-classes/channelpedia_info.csv'
PATH_CHEMBL_ASSAYS_COMPOUND = 'data/chembl-assays/output_compound-assays_v1.1.dat' # Not used
PATH_CHEMBL_ASSAYS_DRUG = 'data/chembl-assays/output_drug-assays_v1.1.dat' # Not used
PATH_CHEMBL_HUMAN_ASSAYS_COMPOUND = 'data/chembl-assays/output_human_compounds-assays_v1.1.dat'
PATH_CHEMBL_HUMAN_ASSAYS_DRUG = 'data/chembl-assays/output_human_drug-assays_v1.1.dat'
PATH_GO_ION_CHANNELS = 'data/channel-classes/ion-channels-1.csv'
PATH_GO_QUICKGO = 'data/go/QuickGO-ion-channel-COMBINED-human.dat'
PATH_GO_AMIGO = 'data/go/go_0006811_taxon_9606.dat'
PATH_HPA_NORMAL = 'data/hpa/normal_tissue.tsv'
PATH_HPA_PATHOLOGY = 'data/hpa/pathology.tsv'
PATH_HPA_RNA_TISSUE = 'data/hpa/rna_tissue.tsv'
PATH_HPA_TRANSLATION_TABLE = 'data/hpa/hpa_tissue_translation_table.csv'
PATH_IN_BETSE = 'data/channel-classes/in_betse.csv'
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
    tissue_name_list = []
    with open(PATH_BTO_TISSUES, 'rU') as tissue_file:
        for line in tissue_file:
            tissue_name_list.append(line.strip())
    bto_parser = go_obo_parser.parseGOOBO(PATH_BTO_BRENDA_TISSUE_OBO)
    for bto_record in bto_parser:
        tissue_name = bto_record['name']
        if tissue_name in tissue_name_list:
            bto_id = bto_record['id']
            db.add_tissue(tissue_name, bto_id, '')


def input_biogps_tissue_data():
    db = icdb.IonChannelDatabase(PATH_DB)
    with open(PATH_BIOGPS_TRANSLATION_TABLE, 'rU') as translation_file:
        translation_reader = csv.reader(translation_file)
        translation_reader.next()
        for row in translation_reader:
            biogps_tissue = row[0]
            bto_tissue = row[1]
            db.add_db_tissue(EXTDB_BIOGPS, bto_tissue, biogps_tissue)


def input_hpa_tissue_data():
    db = icdb.IonChannelDatabase(PATH_DB)
    with open(PATH_HPA_TRANSLATION_TABLE, 'rU') as translation_file:
        translation_reader = csv.reader(translation_file)
        translation_reader.next()
        for row in translation_reader:
            hpa_tissue = row[0]
            bto_tissue = row[1]
            if bto_tissue != '':
                db.add_db_tissue(EXTDB_HPA, bto_tissue, hpa_tissue)


def setup_db_tissue_table():
    db = icdb.IonChannelDatabase(PATH_DB)
    db.create_db_tissue_table()
    input_biogps_tissue_data()
    input_hpa_tissue_data()


# Routine is not being used
def setup_genbank_uniprot_table():
    db = icdb.IonChannelDatabase(PATH_DB)
    db.create_genbank_uniprot_table()
    with open(PATH_UNIPROT_GENBANK, 'rU') as gu_file:
        gu_file.next()
        for line in gu_file:
            line_split = line.split('\t')
            if len(line_split) == 2:
                gbac = line_split[0].strip()
                upac = line_split[1].strip()
                if upac != '' and not db.exists_genbank_uniprot(gbac, upac):
                    db.add_genbank_uniprot(gbac, upac)


def input_chembl_compound_data():
    db = icdb.IonChannelDatabase(PATH_DB)
    with open(PATH_CHEMBL_HUMAN_ASSAYS_COMPOUND, 'rU') as compound_file:
        for line in compound_file:
            line = line.decode('utf_8').encode('utf_8')
            line_split = line.split('\t')
            target_upac = line_split[0].strip()
            target_chembl_id = line_split[1].strip()
            target_type = line_split[2].strip() # Ignored
            target_name = line_split[3].strip()
            compound_chembl_id = line_split[4].strip()
            compound_name = line_split[5].strip()
            iupac_name = line_split[6].strip()
            assay_chembl_id = line_split[7].strip()
            assay_standard_type = line_split[8].strip()
            assay_standard_relation = line_split[9].strip() # Ignored
            assay_value = line_split[10].strip()
            assay_units = line_split[11].strip()
            assay_type = line_split[12].strip()
            if assay_type not in ASSAY_TYPE_CHOICES:
                print 'Unexpected assay_type %s' % (assay_type,)
            assay_description = line_split[13].strip() # Ignored
            if not db.exists_compound_by_chembl_id(compound_chembl_id):
                db.add_compound(
                    compound_name, compound_type=COMPOUND_TYPE_COMPOUND,
                    chembl_id=compound_chembl_id, synonyms=iupac_name,
                    source_db_name=EXTDB_CHEMBL)
            if not db.exists_protein(target_upac):
                db.add_protein(
                    target_upac, name=target_name, chembl_id=target_chembl_id)
            compound_id = db.get_compound_id_by_chembl_id(compound_chembl_id)
            db.add_interaction(
                target_upac, compound_id, assay_value=assay_value,
                assay_units=assay_units,
                assay_standard_type=assay_standard_type,
                assay_type=assay_type, assay_chembl_id=assay_chembl_id,
                source_db_name=EXTDB_CHEMBL)


def input_chembl_drug_data():
    db = icdb.IonChannelDatabase(PATH_DB)
    with open(PATH_CHEMBL_HUMAN_ASSAYS_DRUG, 'rU') as drug_file:
        for line in drug_file:
            line = line.decode('utf_8').encode('utf_8')
            line_split = line.split('\t')
            target_upac = line_split[0].strip()
            target_chembl_id = line_split[1].strip()
            target_type = line_split[2].strip() # Ignored
            target_name = line_split[3].strip()
            drug_mechanism = line_split[4].strip()
            if drug_mechanism not in ACTION_TYPE_CHOICES:
                print 'Unexpected action_type %s' % (drug_mechanism,)
            dosed_compound_chembl_id = line_split[5].strip()
            dosed_compound_name = line_split[6].strip()
            active_compound_chembl_id = line_split[7].strip()
            active_compound_name = line_split[8].strip()
            assay_chembl_id = line_split[9].strip()
            assay_standard_type = line_split[10].strip()
            assay_standard_relation = line_split[11].strip() # Ignored
            assay_value = line_split[12].strip()
            assay_units = line_split[13].strip()
            assay_type = line_split[14].strip()
            assay_description = line_split[15].strip() # Ignored
            if not db.exists_compound_by_chembl_id(active_compound_chembl_id):
                db.add_compound(
                    active_compound_name,  compound_type=COMPOUND_TYPE_DRUG,
                    chembl_id=active_compound_chembl_id,
                    synonyms=dosed_compound_name, source_db_name=EXTDB_CHEMBL)
            if not db.exists_protein(target_upac):
                db.add_protein(
                    target_upac, name=target_name, chembl_id=target_chembl_id)
            compound_id = db.get_compound_id_by_chembl_id(
                active_compound_chembl_id)
            db.add_interaction(
                target_upac, compound_id, action_type=drug_mechanism,
                assay_value=assay_value, assay_units=assay_units,
                assay_standard_type=assay_standard_type, assay_type=assay_type,
                assay_chembl_id=assay_chembl_id, source_db_name=EXTDB_CHEMBL)


def input_chembl_assay_data():
    input_chembl_compound_data()
    input_chembl_drug_data()


def input_ion_channel_data():
    db = icdb.IonChannelDatabase(PATH_DB)
    with open(PATH_GO_ION_CHANNELS, 'rU') as prot_file:
        prot_reader = csv.reader(prot_file)
        prot_reader.next()
        for row in prot_reader:
            row = [item.decode('utf_8').encode('utf_8') for item in row]
            upac = row[0]
            gene_symbol = row[1]
            name = row[2]
            ions = row[3]
            gating = row[4]
            if upac != '':
                if gene_symbol == '-':
                    gene_symbol = ''
                if not db.exists_protein(upac):
                    db.add_protein(
                        upac, gene_symbol=gene_symbol, name=name, ions=ions,
                        gating=gating)
                else:
                    old_protein_record = db.lookup_protein(upac)
                    if old_protein_record is not None:
                        old_gene_symbol = old_protein_record['GeneSymbol']
                        if old_gene_symbol == '':
                            db.update_protein_gene_symbol(upac, gene_symbol)
                        elif old_gene_symbol != gene_symbol:
                            db.update_protein_gene_symbol(upac, gene_symbol)
                            print 'Gene symbol mismatch upac=%s old_gene_symbol=%s gene_symbol=%s' \
                                % (upac, old_gene_symbol, gene_symbol)
                        old_name = old_protein_record['Name']
                        if old_name == '':
                            db.update_protein_name(upac, name)
                        db.update_protein_ions(upac, ions)
                        db.update_protein_gating(upac, gating)
                

def input_in_betse_data():
    db = icdb.IonChannelDatabase(PATH_DB)
    with open(PATH_IN_BETSE, 'rU') as in_betse_file:
        in_betse_reader = csv.reader(in_betse_file)
        in_betse_reader.next()
        for row in in_betse_reader:
            gene_symbol = row[0].strip()
            in_betse = row[1].strip()
            upac_list = db.get_uniprot_accnums_by_gene_symbol(gene_symbol)
            for upac in upac_list:
                db.update_protein_in_betse(upac, in_betse)


def setup_protein_compound_interaction_tables():
    db = icdb.IonChannelDatabase(PATH_DB)
    db.create_protein_table()
    db.create_compound_table()
    db.create_interaction_table()
    input_chembl_assay_data()
    input_ion_channel_data()
    input_in_betse_data()


def get_biogps_probeset_id_to_gene_symbol_dict():
    annot_dict = {}
    with open(PATH_BIOGPS_GNF1H_ANNOT, 'rU') as chip_annot_file:
        chip_annot_reader = csv.reader(chip_annot_file, delimiter='\t')
        chip_annot_reader.next()
        for row in chip_annot_reader:
            probeset_id = row[0]
            gene_symbol = row[6]
            if gene_symbol not in ['', 'obsoleted by Celera']:
                annot_dict[probeset_id] = gene_symbol
        chip_annot_file.close()
        return annot_dict


def input_biogps_expression_data():
    db = icdb.IonChannelDatabase(PATH_DB)
    annot_dict = get_biogps_probeset_id_to_gene_symbol_dict()
    with open(PATH_BIOGPS_GCRMA, 'rU') as microarray_file:
        microarray_reader = csv.reader(microarray_file)
        header = microarray_reader.next()
        tissue_list = header[1:]
        for row in microarray_reader:
            probeset_id = row[0]
            if probeset_id in annot_dict:
                gene_symbol = annot_dict[probeset_id].strip()
                upac_list = db.get_uniprot_accnums_by_gene_symbol(gene_symbol)
                if upac_list:
                    for i in xrange(0, len(tissue_list)):
                        tissue_name = tissue_list[i]
                        expr_level = float(row[i + 1])
                        tissue_bto = db.get_tissue_name_by_db_tissue_name(
                            EXTDB_BIOGPS, tissue_name)
                        if tissue_bto != '':
                            for upac in upac_list:
                                db.add_expression(
                                    tissue_bto, upac, expr_level=expr_level,
                                    assay_type=EXPR_ASSAY_MICROARRAY,
                                    dataset_name=DATASET_BIOGPS,
                                    source_db_name=EXTDB_BIOGPS)


# Function is not used
def get_ensembl_to_uniprot_dict():
    ensembl_to_uniprot = {}
    with open(PATH_UNIPROT_ENSEMBL_IDMAPPING, 'rU') as map_file:
        for line in map_file:
            row = line.strip().split('\t')
            upac = row[0]
            ensg = row[2]
            ensembl_to_uniprot[ensg] = upac
        map_file.close()
        return ensembl_to_uniprot


def input_hpa_rna_expression_data():
    db = icdb.IonChannelDatabase(PATH_DB)
    with open(PATH_HPA_RNA_TISSUE, 'rU') as hpa_file:
        hpa_file.next()
        for line in hpa_file:
            row = line.strip().split('\t')
            if len(row) >= 5 and row[0].strip() != '':
                gene_name = row[1]
                sample = row[2]
                value = float(row[3])
                unit = row[4]
                upac_list = db.get_uniprot_accnums_by_gene_symbol(gene_name)
                for upac in upac_list:
                    tissue_bto = db.get_tissue_name_by_db_tissue_name(
                        EXTDB_HPA, sample)
                    if tissue_bto != '':
                        db.add_expression(
                            tissue_bto, upac, expr_level=value,
                            expr_units=unit, assay_type=EXPR_ASSAY_RNASEQ,
                            dataset_name=DATASET_HPA, source_db_name=EXTDB_HPA)


def input_hpa_cancer_expression_data():
    db = icdb.IonChannelDatabase(PATH_DB)
    with open(PATH_HPA_PATHOLOGY, 'rU') as hpa_file:
        hpa_file.next()
        expr_level_qual_choices = [EXPR_H, EXPR_M, EXPR_L, EXPR_ND]
        for line in hpa_file:
            row = line.strip().split('\t')
            if len(row) >= 7:
                gene_name = row[1]
                cancer = row[2]
                upac_list = db.get_uniprot_accnums_by_gene_symbol(gene_name)
                for upac in upac_list:
                    tissue_bto = db.get_tissue_name_by_db_tissue_name(
                        EXTDB_HPA, cancer)
                    if tissue_bto != '':
                        if not '' in row[3:7]:
                            patient_count_list = map(int, row[3:7])
                            patient_count_max = max(patient_count_list)
                            patient_count_max_index = 0
                            for i in xrange(len(patient_count_list)-1, -1, -1):
                                if patient_count_list[i] == patient_count_max:
                                    patient_count_max_index = i
                            expr_level_qual = expr_level_qual_choices[
                                patient_count_max_index]
                            db.add_expression(
                                tissue_bto, upac,
                                expr_level_qual=expr_level_qual,
                                assay_type=EXPR_ASSAY_IMMUNO,
                                dataset_name=DATASET_HPA,
                                source_db_name=EXTDB_HPA)


def input_hpa_normal_expression_data():
    db = icdb.IonChannelDatabase(PATH_DB)
    with open(PATH_HPA_NORMAL, 'rU') as hpa_file:
        hpa_file.next()
        for line in hpa_file:
            row = line.strip().split('\t')
            if len(row) >= 6:
                gene_name = row[1]
                tissue = row[2]
                level = row[4]
                upac_list = db.get_uniprot_accnums_by_gene_symbol(gene_name)
                for upac in upac_list:
                    tissue_bto = db.get_tissue_name_by_db_tissue_name(
                        EXTDB_HPA, tissue)
                    if tissue_bto != '':
                        if level in [EXPR_H, EXPR_M, EXPR_L, EXPR_ND]:
                            db.add_expression(
                                tissue_bto, upac, expr_level_qual=level,
                                assay_type=EXPR_ASSAY_IMMUNO,
                                dataset_name=DATASET_HPA,
                                source_db_name=EXTDB_HPA)


def setup_expression_table():
    db = icdb.IonChannelDatabase(PATH_DB)
    db.create_expression_table()
    input_biogps_expression_data()
    input_hpa_rna_expression_data()
    input_hpa_cancer_expression_data()
    input_hpa_normal_expression_data()


def get_channelpedia_dict():
    channelpedia_dict = {}
    with open(PATH_CHANNELPEDIA_INFO, 'rU') as channelpedia_file:
        channelpedia_reader = csv.reader(channelpedia_file)
        channelpedia_reader.next()
        for row in channelpedia_reader:
            row = [item.decode('utf_8').encode('utf_8') for item in row]
            subclass = row[1]
            intro_text = row[2]
            url = row[3]
            channelpedia_dict[subclass] = (intro_text, url)
    return channelpedia_dict


def setup_channel_class_tables():
    db = icdb.IonChannelDatabase(PATH_DB)
    db.create_channel_super_class_table()
    db.create_channel_class_table()
    db.create_channel_sub_class_table()
    channelpedia_dict = get_channelpedia_dict()
    with open(PATH_CHANNEL_CLASSES, 'rU') as class_file:
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
                    channelpedia_text, channelpedia_url = channelpedia_dict[
                        channel_subclass]
                else:
                    channelpedia_text = ''
                    channelpedia_url = ''
                db.add_channel_sub_class(
                    channel_subclass, channel_class,
                    channelpedia_text=channelpedia_text,
                    channelpedia_url=channelpedia_url, subfamily=subfamily)
            upac_list = db.get_uniprot_accnums_by_gene_symbol(gene_symbol)
            for upac in upac_list:
                db.update_protein_sub_class(upac, channel_subclass)


#def setup_go_term_table():
#    db = icdb.IonChannelDatabase(PATH_DB)
#    db.create_go_term_table()
#    with open(PATH_GO_QUICKGO, 'rU') as go_file:
#        go_reader = csv.reader(go_file, delimiter='\t')
#        go_reader.next()
#        for row in go_reader:
#            db_name = row[0]
#            upac = row[1]
#            gene_symbol = row[2]
#            go_qualifier = row[3]
#            goid = row[4]
#            go_name = row[5]
#            taxon_id = row[6]
#            gene_product_name = row[7]
#            go_aspect = row[8]
#            if db_name == 'UniProtKB':
#                if db.exists_protein(upac):
#                    old_protein_record = db.lookup_protein(upac)
#                    old_gene_symbol = old_protein_record['GeneSymbol']
#                    if old_gene_symbol == '':
#                        db.update_protein_gene_symbol(upac, gene_symbol)
#                    elif old_gene_symbol != gene_symbol:
#                        db.update_protein_gene_symbol(upac, gene_symbol)
#                    old_name = old_protein_record['Name']
#                    if old_name == '':
#                        db.update_protein_name(upac, gene_product_name)
#                    db.add_go_term(
#                        upac, goid, name=go_name, qualifier=go_qualifier,
#                        aspect=go_aspect)

def setup_go_term_table():
    db = icdb.IonChannelDatabase(PATH_DB)
    db.create_go_term_table()
    with open(PATH_GO_AMIGO, 'rU') as go_file:
        go_reader = csv.reader(go_file, delimiter='\t')
        go_reader.next()
        for row in go_reader:
            bioentity = row[0]
            db_name, upac = bioentity.split(':')
            bioentity_name = row[1]
            qualifier = row[2]
            annotation_class = row[3]
            annotation_extension_json = row[4]
            assigned_by = row[5]
            taxon = row[6]
            evidence_type = row[7]
            evidence_with = row[8]
            panther_family = row[9]
            type = row[10]
            bioentity_isoform = row[11]
            reference = row[12]
            date = row[13]
            if db_name == 'UniProtKB' and db.exists_protein(upac):
                db.add_go_term(upac, 'GO:0006811')


def setup_specificity_table_with_jp_method():
    db = icdb.IonChannelDatabase(PATH_DB)
    db.create_specificity_table()
    upac_list = db.get_protein_uniprot_accnums()
    tissue_name_list = db.get_tissue_names()
    specificity_score_threshold = 5.
    specificity_score_dict = {}
    expr_table = {}
    for upac in upac_list:
        specificity_score = 0.
        tissue_count = 0
        for tissue_name in tissue_name_list:
            expr_level_list = \
                db.get_expression_level_by_uniprot_accnum_tissue_dataset(
                    upac, tissue_name, EXTDB_HPA)
            if expr_level_list:
                expr_level_rna = expr_level_list[0]
                expr_table[(tissue_name, upac)] = expr_level_rna
                if expr_level_rna > specificity_score_threshold:
                    specificity_score += 1.
                tissue_count += 1
        if tissue_count > 0:
            specificity_score = 1. - specificity_score / float(tissue_count)
        else:
            specificity_score = None
        specificity_score_dict[upac] = specificity_score
    for tissue_name, upac in expr_table:
        expr_level_rna = expr_table[(tissue_name, upac)]
        specificity_score_for_tissue = \
            specificity_score_dict[upac] * expr_level_rna
        if specificity_score_for_tissue != None:
            db.add_specificity(tissue_name, upac, specificity_score_for_tissue)


def vacuum_db():
    db = icdb.IonChannelDatabase(PATH_DB)
    db.vacuum_db()


def dump_db():
    db = icdb.IonChannelDatabase(PATH_DB)
    db.dump_db()


def print_db_stats():
    db = icdb.IonChannelDatabase(PATH_DB)
    db.print_db_stats()
