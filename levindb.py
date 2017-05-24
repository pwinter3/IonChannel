import sqlite3
import csv
import go_obo_parser

DB_FILENAME = 'levin.db'
OBO_FILENAME = 'data/bto/BrendaTissueOBO.obo'
PROTEIN_FILENAME = 'data/go/ion-channels.csv'
EXTDB_BIOGPS = 'biogps'
GPL96_FILENAME = 'data/geo/GPL96-57554.txt'
GENBANKUNIPROT_FILENAME = 'data/uniprot/GenBankUniProt.txt' 
BIOGPS_GCRMA = 'data/biogps/U133AGNF1B.gcrma.avg.csv'
BIOGPS_DATASET = 'U133AGNF1B'
BIOGPS_GNF1H_ANNOT = 'data/biogps/gnf1h.annot2007.tsv'
BIOGPS_TRANSLATION_TABLE = 'data/biogps/biogps_translation_table.csv'
DRUGBANK_FILENAME = 'data/target-compound/output-organized-Drugbank.dat'
CHEMBL_COMPOUND_FILENAME = \
        'data/target-compound/output-query_ChEMBL-uniprot-compound.dat'

def create_tables():
    conn = sqlite3.connect(DB_FILENAME)
    curs = conn.cursor()
    curs.execute('''CREATE TABLE Tissue(
        Name TEXT PRIMARY KEY NOT NULL,
        BTOId TEXT UNIQUE NOT NULL,
        ParentTissueName TEXT NOT NULL
        );''')
    curs.execute('''CREATE TABLE Protein(
        UniProtAccNum TEXT PRIMARY KEY NOT NULL,
        GeneSymbol TEXT UNIQUE NOT NULL,
        Name TEXT NOT NULL,
        ProcessFunction TEXT NOT NULL
        );''')
    curs.execute('''CREATE TABLE ExternalDB(
        Name TEXT PRIMARY KEY NOT NULL,
        URL TEXT NOT NULL
        );''')
    curs.execute('''CREATE TABLE DBTissue(
        ExternalDBName TEXT NOT NULL,
        TissueName TEXT NOT NULL,
        DBEquivalentTissueName TEXT NOT NULL,
        PRIMARY KEY (ExternalDBName, DBEquivalentTissueName)
        );''')
    curs.execute('''CREATE TABLE Expression(
        Id INTEGER PRIMARY KEY AUTOINCREMENT,
        TissueName TEXT NOT NULL,
        ProteinUniProtAccNum TEXT NOT NULL,
        ExprLevel REAL NOT NULL,
        ExprUnits TEXT NOT NULL,
        AssayType TEXT NOT NULL,
        DatasetName TEXT NOT NULL,
        SourceDBName TEXT NOT NULL
        );''')
    curs.execute('''CREATE TABLE DBGene(
        ExternalDBName TEXT NOT NULL,
        GenbankAccNum TEXT NOT NULL,
        ProbeID TEXT NOT NULL,
        PRIMARY KEY (ExternalDBName, GenbankAccNum)
        );''')
    curs.execute('''CREATE TABLE GenbankUniprot(
        GenbankAccNum TEXT NOT NULL,
        UniProtAccNum TEXT NOT NULL,
        PRIMARY KEY (GenbankAccNum, UniprotAccNum)
        );''')
    curs.execute('''CREATE TABLE Compound(
        Id INTEGER PRIMARY KEY AUTOINCREMENT,
        SMILES TEXT NOT NULL,
        InChI TEXT NOT NULL,
        Name TEXT NOT NULL,
        Synonyms TEXT NOT NULL,
        ApprovalStatus TEXT NOT NULL,
        FirstApprovalYear TEXT NOT NULL,
        SourceDBName NOT NULL
        );''')
    curs.execute('''CREATE TABLE Interaction(
        Id INTEGER PRIMARY KEY AUTOINCREMENT,
        TargetUniProtAccNum TEXT NOT NULL,
        CompoundId INTEGER NOT NULL,
        ActionType TEXT NOT NULL,
        ActionDesc TEXT NOT NULL,
        Strength REAL NOT NULL,
        StrengthUnits TEXT NOT NULL,
        AssayType TEXT NOT NULL,
        SourceDBName TEXT NOT NULL
        );''')
    conn.commit()
    conn.close()

def add_externaldb(name, url):
    conn = sqlite3.connect(DB_FILENAME)
    curs = conn.cursor()
    curs.execute('''INSERT INTO ExternalDB (Name, URL)
        VALUES (?, ?)''', (name, url))
    conn.commit()
    conn.close()

def add_tissue(tissue_name, bto_id='', parent_tissue_name=''):
    conn = sqlite3.connect(DB_FILENAME)
    curs = conn.cursor()
    curs.execute('''INSERT INTO Tissue (Name, BTOId, ParentTissueName)
        VALUES (?, ?, ?)''', (tissue_name, bto_id, parent_tissue_name))
    conn.commit()
    conn.close()

def add_protein(uniprot_accnum, gene_symbol='', name='', process_function=''):
    conn = sqlite3.connect(DB_FILENAME)
    curs = conn.cursor()
    curs.execute('''INSERT INTO Protein (UniProtAccNum, GeneSymbol,
        Name, ProcessFunction)
        VALUES (?, ?, ?, ?)''',
        (uniprot_accnum, gene_symbol, name, process_function))
    conn.commit()
    conn.close()

def add_expression(tissue_name, uniprot_accnum, expr_level, expr_units='',
        assay_type='', dataset_name='', sourcedb_name=''):
    conn = sqlite3.connect(DB_FILENAME)
    curs = conn.cursor()
    curs.execute('''INSERT INTO Expression (TissueName, ProteinUniProtAccNum,
        ExprLevel, ExprUnits, AssayType, DatasetName, SourceDBName)
        VALUES (?, ?, ?, ?, ?, ?, ?)''', (tissue_name, uniprot_accnum,
        expr_level, expr_units, assay_type, dataset_name, sourcedb_name))
    conn.commit()
    conn.close()

def add_dbgene(externaldb_name, genbank_accnum, probe_id):
    conn = sqlite3.connect(DB_FILENAME)
    curs = conn.cursor()
    curs.execute('''INSERT INTO DBGene (ExternalDBName, GenbankAccNum, ProbeID)
        VALUES (?, ?, ?)''', (externaldb_name, genbank_accnum, probe_id))
    conn.commit()
    conn.close()

def add_genbankuniprot(genbank_accnum, uniprot_accnum):
    conn = sqlite3.connect(DB_FILENAME)
    curs = conn.cursor()
    curs.execute('''INSERT INTO GenbankUniprot (GenbankAccNum, UniProtAccNum)
        VALUES (?, ?)''', (genbank_accnum, uniprot_accnum))
    conn.commit()
    conn.close()    

def add_compound(smiles, inchi, name, synonyms='', approval_status='',
        first_approval_year='', sourcedb_name=''):
    conn = sqlite3.connect(DB_FILENAME)
    curs = conn.cursor()
    curs.execute('''INSERT INTO Compound (SMILES, InChI, Name, Synonyms,
        ApprovalStatus, FirstApprovalYear, SourceDBName)
        VALUES (?, ?, ?, ?, ?, ?, ?)''',
        (smiles, inchi, name, synonyms, approval_status, first_approval_year, 
        sourcedb_name))
    conn.commit()
    conn.close()    

def add_interaction(target_uniprot_accnum, compound_id, action_type='',
        action_desc='', strength=0, strength_units='', assay_type='',
        sourcedb_name = ''):
    conn = sqlite3.connect(DB_FILENAME)
    curs = conn.cursor()
    curs.execute('''INSERT INTO Interaction (TargetUniProtAccNum, CompoundId,
        ActionType, ActionDesc, Strength, StrengthUnits, AssayType, 
        SourceDBName)
        VALUES (?, ?, ?, ?, ?, ?, ?, ?)''',
        (target_uniprot_accnum, compound_id, action_type, action_desc, 
        strength, strength_units, assay_type, sourcedb_name))
    conn.commit()
    conn.close()    

def add_dbtissue(externaldb_name, tissue_name, db_equivalent_tissue_name):
    conn = sqlite3.connect(DB_FILENAME)
    curs = conn.cursor()
    curs.execute('''INSERT INTO DBTissue (ExternalDBName, TissueName,
        DBEquivalentTissueName)
        VALUES (?, ?, ?)''', (externaldb_name, tissue_name, 
        db_equivalent_tissue_name))
    conn.commit()
    conn.close()

def exists_genbankuniprot(genbank_accnum, uniprot_accnum):
    conn = sqlite3.connect(DB_FILENAME)
    curs = conn.cursor()
    expr_list = []
    curs.execute('''SELECT EXISTS(SELECT 1 FROM GenbankUniprot
            WHERE (GenbankAccNum=? AND UniProtAccNum=?) LIMIT 1)''',
            (genbank_accnum, uniprot_accnum))
    row = curs.fetchone()
    if row[0] == 1:
        return True
    else:
        return False

def exists_protein(uniprot_accnum):
    conn = sqlite3.connect(DB_FILENAME)
    curs = conn.cursor()
    expr_list = []
    curs.execute('''SELECT EXISTS(SELECT 1 FROM Protein
            WHERE UniProtAccNum=? LIMIT 1)''',
            (uniprot_accnum,))
    row = curs.fetchone()
    if row[0] == 1:
        return True
    else:
        return False

def exists_compound(compound_name):
    conn = sqlite3.connect(DB_FILENAME)
    curs = conn.cursor()
    expr_list = []
    curs.execute('''SELECT EXISTS(SELECT 1 FROM Compound
            WHERE Name=? LIMIT 1)''',
            (compound_name,))
    row = curs.fetchone()
    if row[0] == 1:
        return True
    else:
        return False

def exists_interaction(uniprot_accnum, compound_id):
    conn = sqlite3.connect(DB_FILENAME)
    curs = conn.cursor()
    expr_list = []
    curs.execute('''SELECT EXISTS(SELECT 1 FROM Interaction
            WHERE TargetUniProtAccNum=? AND CompoundId=? LIMIT 1)''',
            (uniprot_accnum, compound_id))
    row = curs.fetchone()
    if row[0] == 1:
        return True
    else:
        return False

def dump_proteins():
    conn = sqlite3.connect(DB_FILENAME)
    curs = conn.cursor()
    for row in curs.execute('''SELECT * FROM Protein'''):
        print row

def dump_tissues():
    conn = sqlite3.connect(DB_FILENAME)
    curs = conn.cursor()
    for row in curs.execute('''SELECT * FROM Tissue'''):
        print row

def dump_database():
    conn = sqlite3.connect(DB_FILENAME)
    with open('dump.sql', 'w') as f:
        for line in conn.iterdump():
            f.write('%s\n' % line)

def lookup_uniprot_accnum_by_gene_symbol(gene_symbol):
    conn = sqlite3.connect(DB_FILENAME)
    curs = conn.cursor()
    curs.execute('''SELECT UniProtAccNum FROM Protein WHERE GeneSymbol=?''', 
        (gene_symbol,))
    resultset = curs.fetchone()
    if resultset == None:
        result = ''
    else:
        result = resultset[0]
    conn.close()
    return result

def lookup_gene_symbol_by_uniprot_accnum(uniprot_accnum):
    conn = sqlite3.connect(DB_FILENAME)
    curs = conn.cursor()
    curs.execute('''SELECT GeneSymbol FROM Protein WHERE UniProtAccNum=?''', 
        (uniprot_accnum,))
    resultset = curs.fetchone()
    if resultset == None:
        result = ''
    else:
        result = resultset[0]
    conn.close()
    return result

def lookup_dbtissue(externaldb_name, db_equivalent_tissue_name):
    conn = sqlite3.connect(DB_FILENAME)
    curs = conn.cursor()
    curs.execute('''SELECT TissueName FROM DBTissue WHERE ExternalDBName=?
        AND DBEquivalentTissueName=?''', 
        (externaldb_name, db_equivalent_tissue_name))
    resultset = curs.fetchone()
    if resultset == None:
        result = ''
    else:
        result = resultset[0]
    conn.close()
    return result

def lookup_compound_by_name(compound_name):
    conn = sqlite3.connect(DB_FILENAME)
    curs = conn.cursor()
    curs.execute('''SELECT Id FROM Compound WHERE Name=?''', 
        (compound_name,))
    resultset = curs.fetchone()
    if resultset == None:
        result = ''
    else:
        result = resultset[0]
    conn.close()
    return result

def setup_externaldb():
    add_externaldb('biogps', 'http://biogps.org')
    add_externaldb('chembl', 'https://www.ebi.ac.uk/chembl/')
    add_externaldb('brenda', 'http://www.brenda-enzymes.org')

def cleanup_externaldb():
    conn = sqlite3.connect(DB_FILENAME)
    curs = conn.cursor()
    curs.execute('''DELETE FROM ExternalDB''')
    conn.commit()
    conn.close()

def setup_ontology():
    parser = go_obo_parser.parseGOOBO(OBO_FILENAME)
    for record in parser:
        bto_id = record['id']
        name = record['name']
        add_tissue(name, bto_id=bto_id)

def cleanup_ontology():
    conn = sqlite3.connect(DB_FILENAME)
    curs = conn.cursor()
    curs.execute('''DELETE FROM Tissue''')
    conn.commit()
    conn.close()

def setup_biogps():
    add_db('biogps', 'http://biogps.org/')
    translation_f = open(BIOGPS_TRANSLATION_TABLE, 'rU')
    translation_r = csv.reader(translation_f)
    translation_r.next()
    for row in translation_r:
        biogps_tissue = row[0]
        bto_tissue = row[1]
        add_dbtissue('biogps', bto_tissue, biogps_tissue)
    translation_f.close()

def cleanup_dbtissue():
    conn = sqlite3.connect(DB_FILENAME)
    curs = conn.cursor()
    curs.execute('''DELETE FROM DBTissue''')
    conn.commit()
    conn.close()

def setup_proteins():
    prot_file = open(PROTEIN_FILENAME, 'rU')
    prot_reader = csv.reader(prot_file)
    prot_reader.next()
    for row in prot_reader:
        name = row[6]
        process_function = row[7]
        #if process_function.startswith('NOT AN ION CHANNEL'):
        #    continue
        uniprot_accnum = row[4]
        if uniprot_accnum == '':
            continue
        if not uniprot_accnum[0] in ['O', 'P', 'Q']:
            continue
        gene_symbol = row[5]
        if gene_symbol == '-':
            gene_symbol = ''
        if gene_symbol == '':
            continue
        check_uniprot_accnum = lookup_uniprot_accnum_by_gene_symbol(gene_symbol)
        if check_uniprot_accnum != '':
            continue
        check_gene_symbol = lookup_gene_symbol_by_uniprot_accnum(uniprot_accnum)
        if check_gene_symbol != '':
            continue
        add_protein(uniprot_accnum, gene_symbol, name, process_function)
    prot_file.close()

def cleanup_proteins():
    conn = sqlite3.connect(DB_FILENAME)
    curs = conn.cursor()
    curs.execute('''DELETE FROM Protein''')
    conn.commit()
    conn.close()

def setup_genbankuniprot():
    gu_file = open(GENBANKUNIPROT_FILENAME, 'rU')
    line_num = 0
    for line in gu_file:
        if line_num != 0:
            line_split = line.split('\t')
            if len(line_split) == 2:
                genbank_accnum = line_split[0].strip()
                uniprot_accnum = line_split[1].strip()
                if uniprot_accnum != '':
                    if not exists_genbankuniprot(genbank_accnum, 
                            uniprot_accnum):
                        try:
                            add_genbankuniprot(genbank_accnum, uniprot_accnum)
                        except Exception as e:
                            print e
                            print genbank_accnum, uniprot_accnum
        line_num += 1
    gu_file.close()

def cleanup_genbankuniprot():
    conn = sqlite3.connect(DB_FILENAME)
    curs = conn.cursor()
    curs.execute('''DELETE FROM GenbankUniprot''')
    conn.commit()
    conn.close()

def setup_chembl():
    chembl_compound_file = open(CHEMBL_COMPOUND_FILENAME)
    chembl_compound_file.next()
    for line in chembl_compound_file:
        line = ''.join([i if ord(i) < 128 else '' for i in line])
        uniprot = line[0:14].strip()
        #target_name = line[14:96].strip()
        #target_id = line[96:114].strip()
        #compound_key = line[114:128].strip()
        chembl_compound_id = line[128:146].strip()
        compound_name_in_doc = line[146:348].strip()
        first_approval = line[348:356].strip()
        withdrawn_flag = line[356:].strip()
        # This is a temporary way of checking if the compound exists
        if exists_protein(uniprot) and not exists_compound(chembl_compound_id):
            if withdrawn_flag == 1:
                approval_status = 'WITHDRAWN'
            else:
                approval_status = 'UNKNOWN'
            add_compound('', '', chembl_compound_id, compound_name_in_doc,
                    approval_status, first_approval, 'chembl')
            print 'Added compound %s' % chembl_compound_id
        if exists_protein(uniprot) and exists_compound(chembl_compound_id):
            compound_id = lookup_compound_by_name(chembl_compound_id)
            if not exists_interaction(uniprot, compound_id):
                add_interaction(uniprot, compound_id, sourcedb_name='chembl')
                print 'Added interaction %s %s' % (uniprot, chembl_compound_id)
    chembl_compound_file.close()

def setup_hmdb():
    pass

def setup_drugbank():
    pass

def cleanup_compounds():
    conn = sqlite3.connect(DB_FILENAME)
    curs = conn.cursor()
    curs.execute('''DELETE FROM Compound''')
    conn.commit()
    conn.close()

def get_expression(tissue_name, uniprot_accnum):
    conn = sqlite3.connect(DB_FILENAME)
    curs = conn.cursor()
    expr_list = []
    for row in curs.execute('''SELECT ExprLevel FROM Expression 
            WHERE TissueName = ? AND ProteinUniProtAccNum = ?''',
            (tissue_name, uniprot_accnum)):
        print row
        expr_list.append(row[0])

def get_expression_values(tissue_list, protein_list):
    conn = sqlite3.connect(DB_FILENAME)
    curs = conn.cursor() 
    full_list = []
    for row in curs.execute('''SELECT * FROM Expression 
            WHERE TissueName = ? AND ProteinUniProtAccNum = ?''',
            (tissue_list, protein_list)):
        print row
        full_list.append(row[0])

def process_biogps_annotation():
    chip_annot_f = open(BIOGPS_GNF1H_ANNOT, 'rU')
    chip_annot_r = csv.reader(chip_annot_f, delimiter='\t')
    chip_annot_r.next()
    annot_dict = {}
    for row in chip_annot_r:
        probeset_id = row[0]
        symbol = row[6]
        #if probeset_id.startswith('gnf1h'):
        #    probeset_id = probeset_id[5:].lstrip('0')
        if symbol != '':
            annot_dict[probeset_id] = symbol
    chip_annot_f.close()
    return annot_dict

def process_biogps_microarray(annot_dict):
    microarray_f = open(BIOGPS_GCRMA, 'rU')
    microarray_r = csv.reader(microarray_f)
    header = microarray_r.next()
    tissue_list = header[1:]
    for row in microarray_r:
        probeset_id = row[0]
        if not probeset_id in annot_dict:
            continue
        symbol = annot_dict[probeset_id]
        uniprot_accnum = lookup_uniprot_accnum_by_gene_symbol(symbol.strip())
        if uniprot_accnum == '':
            continue
        for i in xrange(0, len(tissue_list)):
            expr_level = float(row[i + 1])
            tissue_name = lookup_dbtissue('biogps', tissue_list[i])
            print tissue_name
            add_expression(tissue_name, uniprot_accnum, expr_level, 
                    assay_type='microarray', dataset_name=BIOGPS_DATASET, 
                    sourcedb_name='biogps')
    microarray_f.close()

def process_biogps_data():
    annot_dict = process_biogps_annotation()
    process_biogps_microarray(annot_dict)

def cleanup_expression():
    conn = sqlite3.connect(DB_FILENAME)
    curs = conn.cursor()
    curs.execute('''DELETE FROM Expression''')
    conn.commit()
    conn.close()


