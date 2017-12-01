import sqlite3
import csv
import unicodecsv
import go_obo_parser
import os
import codecs

DB_FILENAME = 'levin.db'
OBO_FILENAME = 'data/bto/BrendaTissueOBO.obo'
TISSUE_FILENAME = 'data/bto/tissues.txt'
PROTEIN_FILENAME = 'data/go/ion-channels-1.csv'
GPL96_FILENAME = 'data/geo/GPL96-57554.txt'
GENBANKUNIPROT_FILENAME = 'data/uniprot/GenBankUniProt.txt' 
BIOGPS_GCRMA = 'data/biogps/U133AGNF1B.gcrma.avg.csv'
BIOGPS_GNF1H_ANNOT = 'data/biogps/gnf1h.annot2007.tsv'
BIOGPS_TRANSLATION_TABLE = 'data/biogps/biogps_translation_table.csv'
DRUGBANK_FILENAME = 'data/target-compound/output-organized-Drugbank.dat'
CHEMBL_COMPOUND_FILENAME = \
        'data/target-compound/output-query_ChEMBL-uniprot-compound.dat'
TARGET_COMPOUND_FILENAME = 'data/target-compound/target-compound.csv'
CHANNEL_CLASSES_FILENAME = 'data/channel-classes/channel-classes.csv'
HUMAN_IDMAPPING = 'data/uniprot/HUMAN_9606_idmapping_Gene_Name.dat'

OUTPUT_ORGANIZED_DRUG_INFO_DRUGBANK = 'data/Target-CompoundTables-2017-01/output-organized_drug-info-Drugbank.dat'
OUTPUT_ORGANIZED_DRUG_INFO_HMDB = 'data/Target-CompoundTables-2017-01/output-organized_drug-info-HMDB.dat'
OUTPUT_ORGANIZED_TARGET_INFO_DRUGBANK = 'data/Target-CompoundTables-2017-01/output-organized_target-info-Drugbank.dat'
OUTPUT_ORGANIZED_TARGET_INFO_HMDB = 'data/Target-CompoundTables-2017-01/output-organized_target-info-HMDB.dat'
OUTPUT_ORGANIZED_TTD = 'data/Target-CompoundTables-2017-01/output-organized-TTD.dat'
OUTPUT_QUERY_CHEMBL_UNIPROT_COMPOUND = 'data/Target-CompoundTables-2017-01/output-query_ChEMBL-uniprot-compound.dat'
OUTPUT_QUERY_CHEMBL_UNIPROT_DRUG = 'data/Target-CompoundTables-2017-01/output-query_ChEMBL-uniprot-drug.dat'
OUTPUT_COMPOUND_ASSAYS = 'data/Chembl-query_2017-08-16/output_compound-assays_v1.1.dat'
OUTPUT_DRUG_ASSAYS = 'data/Chembl-query_2017-08-16/output_drug-assays_v1.1.dat'

HPA_RNA_TISSUE = 'data/hpa/rna_tissue.tsv'
UNIPROT_ENSEMBL_MAPPING = 'data/uniprot/HUMAN_9606_idmapping_Ensembl.dat'
HPA_TISSUE_TRANSLATION_TABLE = 'data/hpa/hpa_tissue_translation_table.csv'

EXTDB_BIOGPS = 'biogps'
EXTDB_CHEMBL = 'chembl'
EXTDB_HPA = 'hpa'
BIOGPS_DATASET = 'U133AGNF1B'
HPA_DATASET = 'HPA_rna_tissue'

def create_tables():
    conn = sqlite3.connect(DB_FILENAME)
    curs = conn.cursor()
    curs.execute('''CREATE TABLE ChannelSuperClass(
        Name TEXT PRIMARY KEY NOT NULL
        );''')
    curs.execute('''CREATE TABLE ChannelClass(
        Name TEXT PRIMARY KEY NOT NULL,
        SuperClass TEXT NOT NULL
        );''')
    curs.execute('''CREATE TABLE ChannelSubClass(
        Name TEXT PRIMARY KEY NOT NULL,
        Class TEXT NOT NULL
        );''')
    curs.execute('''CREATE TABLE Tissue(
        Name TEXT PRIMARY KEY NOT NULL,
        BTOId TEXT UNIQUE NOT NULL,
        ParentTissueName TEXT NOT NULL
        );''')
    curs.execute('''CREATE TABLE Protein(
        UniProtAccNum TEXT PRIMARY KEY NOT NULL,
        GeneSymbol TEXT NOT NULL,
        Name TEXT NOT NULL,
        ProcessFunction TEXT NOT NULL,
        Ions TEXT NOT NULL,
        Gating TEXT NOT NULL,
        InBETSE TEXT NOT NULL,
        IonChannelClassDesc TEXT NOT NULL,
        IonChannelSubClass TEXT NOT NULL,
        ChemblId TEXT NOT NULL
        );''')
    curs.execute('''CREATE TABLE PDB(
        Id INTEGER PRIMARY KEY AUTOINCREMENT NOT NULL,
        UniProtAccNum TEXT NOT NULL,
        PDBID TEXT NOT NULL,
        HumanOnly TEXT NOT NULL
        );''')
    curs.execute('''CREATE TABLE ExternalDB(
        Name TEXT PRIMARY KEY NOT NULL,
        URL TEXT NOT NULL
        );''')
    curs.execute('''CREATE TABLE DBTissue(
        Id INTEGER PRIMARY KEY AUTOINCREMENT NOT NULL,
        ExternalDBName TEXT NOT NULL,
        TissueName TEXT NOT NULL,
        DBEquivalentTissueName TEXT NOT NULL
        );''')
    curs.execute('''CREATE TABLE Expression(
        Id INTEGER PRIMARY KEY AUTOINCREMENT NOT NULL,
        TissueName TEXT NOT NULL,
        ProteinUniProtAccNum TEXT NOT NULL,
        ExprLevel REAL NOT NULL,
        ExprUnits TEXT NOT NULL,
        AssayType TEXT NOT NULL,
        DatasetName TEXT NOT NULL,
        SourceDBName TEXT NOT NULL
        );''')
    curs.execute('''CREATE TABLE DBGene(
        Id INTEGER PRIMARY KEY AUTOINCREMENT NOT NULL,
        ExternalDBName TEXT NOT NULL,
        GenbankAccNum TEXT NOT NULL,
        ProbeID TEXT NOT NULL
        );''')
    curs.execute('''CREATE TABLE GenbankUniprot(
        Id INTEGER PRIMARY KEY AUTOINCREMENT NOT NULL,
        GenbankAccNum TEXT NOT NULL,
        UniProtAccNum TEXT NOT NULL
        );''')
    curs.execute('''CREATE TABLE Compound(
        Id INTEGER PRIMARY KEY AUTOINCREMENT NOT NULL,
        SMILES TEXT NOT NULL,
        InChI TEXT NOT NULL,
        Name TEXT NOT NULL,
        ChemblId TEXT NOT NULL,
        Synonyms TEXT NOT NULL,
        ApprovalStatus TEXT NOT NULL,
        FirstApprovalYear TEXT NOT NULL,
        SourceDBName NOT NULL
        );''')
    curs.execute('''CREATE TABLE Interaction(
        Id INTEGER PRIMARY KEY AUTOINCREMENT NOT NULL,
        TargetUniProtAccNum TEXT NOT NULL,
        CompoundId INTEGER NOT NULL,
        ActionType TEXT NOT NULL,
        ActionDesc TEXT NOT NULL,
        Strength REAL NOT NULL,
        StrengthUnits TEXT NOT NULL,
        AssayType TEXT NOT NULL,
        ChemblId TEXT NOT NULL,
        SourceDBName TEXT NOT NULL
        );''')
    conn.commit()
    conn.close()

def add_channel_superclass(name):
    conn = sqlite3.connect(DB_FILENAME)
    curs = conn.cursor()
    curs.execute('''INSERT INTO ChannelSuperClass (Name)
        VALUES (?)''', (name,))
    conn.commit()
    conn.close()

def add_channel_class(name, superclass):
    conn = sqlite3.connect(DB_FILENAME)
    curs = conn.cursor()
    curs.execute('''INSERT INTO ChannelClass (Name, SuperClass)
        VALUES (?, ?)''', (name, superclass))
    conn.commit()
    conn.close()

def add_channel_subclass(name, channel_class):
    conn = sqlite3.connect(DB_FILENAME)
    curs = conn.cursor()
    curs.execute('''INSERT INTO ChannelSubClass (Name, Class)
        VALUES (?, ?)''', (name, channel_class))
    conn.commit()
    conn.close()

def update_protein_subclass(gene_symbol, subclass):
    conn = sqlite3.connect(DB_FILENAME)
    curs = conn.cursor()
    curs.execute('''UPDATE Protein SET IonChannelSubClass = ?
        WHERE GeneSymbol = ?''', (subclass, gene_symbol))
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

def add_protein(uniprot_accnum, gene_symbol='', name='', process_function='',
        ions='', gating='', in_betse='', ion_channel_class_desc='',
        ion_channel_subclass='', chembl_id=''):
    conn = sqlite3.connect(DB_FILENAME)
    curs = conn.cursor()
    curs.execute('''INSERT INTO Protein (UniProtAccNum, GeneSymbol,
        Name, ProcessFunction, Ions, Gating, InBETSE, IonChannelClassDesc,
        IonChannelSubClass, ChemblId)
        VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?)''',
        (uniprot_accnum, gene_symbol, name, process_function, ions, gating,
        in_betse, ion_channel_class_desc, ion_channel_subclass, chembl_id))
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

def add_compound(smiles, inchi, name, chembl_id='', synonyms='',
        approval_status='', first_approval_year='', sourcedb_name=''):
    conn = sqlite3.connect(DB_FILENAME)
    curs = conn.cursor()
    curs.execute('''INSERT INTO Compound (SMILES, InChI, Name, ChemblId,
        Synonyms, ApprovalStatus, FirstApprovalYear, SourceDBName)
        VALUES (?, ?, ?, ?, ?, ?, ?, ?)''',
        (smiles, inchi, name, chembl_id, synonyms, approval_status,
        first_approval_year, sourcedb_name))
    conn.commit()
    conn.close()    

def add_interaction(target_uniprot_accnum, compound_id, action_type='',
        action_desc='', strength=0, strength_units='', assay_type='',
        chembl_id = '', sourcedb_name = ''):
    conn = sqlite3.connect(DB_FILENAME)
    curs = conn.cursor()
    curs.execute('''INSERT INTO Interaction (TargetUniProtAccNum, CompoundId,
        ActionType, ActionDesc, Strength, StrengthUnits, AssayType, 
        ChemblId, SourceDBName)
        VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?)''',
        (target_uniprot_accnum, compound_id, action_type, action_desc, 
        strength, strength_units, assay_type, chembl_id, sourcedb_name))
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

def exists_channel_superclass(name):
    conn = sqlite3.connect(DB_FILENAME)
    curs = conn.cursor()
    expr_list = []
    curs.execute('''SELECT EXISTS(SELECT 1 FROM ChannelSuperClass
            WHERE Name=? LIMIT 1)''',
            (name,))
    row = curs.fetchone()
    if row[0] == 1:
        return True
    else:
        return False

def exists_channel_class(name):
    conn = sqlite3.connect(DB_FILENAME)
    curs = conn.cursor()
    expr_list = []
    curs.execute('''SELECT EXISTS(SELECT 1 FROM ChannelClass
            WHERE Name=? LIMIT 1)''',
            (name,))
    row = curs.fetchone()
    if row[0] == 1:
        return True
    else:
        return False

def exists_channel_subclass(name):
    conn = sqlite3.connect(DB_FILENAME)
    curs = conn.cursor()
    expr_list = []
    curs.execute('''SELECT EXISTS(SELECT 1 FROM ChannelSubClass
            WHERE Name=? LIMIT 1)''',
            (name,))
    row = curs.fetchone()
    if row[0] == 1:
        return True
    else:
        return False

def exists_compound(chembl_id):
    conn = sqlite3.connect(DB_FILENAME)
    curs = conn.cursor()
    expr_list = []
    curs.execute('''SELECT EXISTS(SELECT 1 FROM Compound
            WHERE ChemblId=? LIMIT 1)''',
            (chembl_id,))
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

def find_prots_where_expr_and_inter():
    conn = sqlite3.connect(DB_FILENAME)
    curs = conn.cursor()
    curs.execute('''SELECT WHERE (EXISTS(SELECT 1 FROM Interaction
            WHERE TargetUniProtAccNum=? AND CompoundId=? LIMIT 1) AND 
            (EXISTS(SELECT 1 FROM Expression
            WHERE ProteinUniprotAccNum=? AND CompoundId=? LIMIT 1))''')

def exists_any_interaction(uniprot_accnum):
    conn = sqlite3.connect(DB_FILENAME)
    curs = conn.cursor()
    curs.execute('''SELECT EXISTS(SELECT 1 FROM Interaction
            WHERE TargetUniProtAccNum=? LIMIT 1)''',
            (uniprot_accnum,))
    row = curs.fetchone()
    if row[0] == 1:
        return True
    else:
        return False

def exists_expression(upac):
    conn = sqlite3.connect(DB_FILENAME)
    curs = conn.cursor()
    curs.execute('''SELECT EXISTS(SELECT 1 FROM Expression
            WHERE ProteinUniProtAccNum=? LIMIT 1)''',
            (upac,))
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
    with codecs.open('dump.sql', 'w', 'utf-8') as f:
        for line in conn.iterdump():
            f.write('%s\n' % line)

def lookup_protein_by_accnum(accnum):
    conn = sqlite3.connect(DB_FILENAME)
    curs = conn.cursor()
    curs.execute('''SELECT UniProtAccNum, GeneSymbol, Name, ProcessFunction,
            Ions, Gating, InBETSE, IonChannelClassDesc, IonChannelSubClass,
            ChemblId FROM Protein WHERE UniProtAccNum=?''', 
        (accnum,))
    resultset = curs.fetchone()
    if resultset == None:
        result = ''
    else:
        result = resultset
    conn.close()
    return result

def is_in_betse(accnum):
    conn = sqlite3.connect(DB_FILENAME)
    curs = conn.cursor()
    curs.execute('''SELECT InBETSE FROM Protein WHERE UniProtAccNum=?''', 
        (accnum,))
    resultset = curs.fetchone()
    if resultset == None:
        result = 'N'
    else:
        result = resultset[0]
    conn.close()
    return result

def lookup_interaction_ids_by_uniprot(accnum):
    conn = sqlite3.connect(DB_FILENAME)
    curs = conn.cursor()
    curs.execute('''SELECT Id FROM Interaction WHERE TargetUniProtAccNum=?''', 
        (accnum,))
    resultset = curs.fetchall()
    conn.close()
    return resultset

def lookup_expression_ids_by_uniprot(accnum):
    conn = sqlite3.connect(DB_FILENAME)
    curs = conn.cursor()
    curs.execute('''SELECT Id FROM Expression WHERE ProteinUniProtAccNum=?''', 
        (accnum,))
    resultset = curs.fetchall()
    conn.close()
    return resultset

def lookup_expression(upac, tissue):
    conn = sqlite3.connect(DB_FILENAME)
    curs = conn.cursor()
    curs.execute('''SELECT ExprLevel FROM Expression WHERE ProteinUniProtAccNum=? AND TissueName=?''', 
        (upac, tissue))
    resultset = curs.fetchall()
    conn.close()
    return resultset

def get_all_protein_uniprot():
    conn = sqlite3.connect(DB_FILENAME)
    curs = conn.cursor()
    curs.execute('''SELECT UniProtAccNum FROM Protein''')
    resultset = curs.fetchall()
    conn.close()
    return resultset

def get_all_tissues():
    conn = sqlite3.connect(DB_FILENAME)
    curs = conn.cursor()
    curs.execute('''SELECT Name FROM Tissue''')
    resultset = curs.fetchall()
    conn.close()
    return resultset

def get_all_expression():
    conn = sqlite3.connect(DB_FILENAME)
    curs = conn.cursor()
    curs.execute('''SELECT Id, TissueName, ProteinUniprotAccNum, ExprLevel FROM Expression''')
    resultset = curs.fetchall()
    conn.close()
    return resultset
    
def get_interaction_by_uniprot(upac):
    conn = sqlite3.connect(DB_FILENAME)
    curs = conn.cursor()
    curs.execute('''SELECT Id FROM Interaction WHERE TargetUniprotAccNum = ?''', (upac,))
    resultset = curs.fetchall()
    conn.close()
    return resultset

def update_protein_name(uniprot_accnum, name):
    conn = sqlite3.connect(DB_FILENAME)
    curs = conn.cursor()
    curs.execute('''UPDATE Protein SET Name = ? WHERE UniProtAccNum = ?''',
        (name, uniprot_accnum))
    conn.commit()
    conn.close()

def update_protein_chembl_id(uniprot_accnum, chembl_id):
    conn = sqlite3.connect(DB_FILENAME)
    curs = conn.cursor()
    curs.execute('''UPDATE Protein SET ChemblId = ? WHERE UniProtAccNum = ?''',
        (chembl_id, uniprot_accnum))
    conn.commit()
    conn.close()

def update_protein_gene_symbol(uniprot_accnum, gene_symbol):
    conn = sqlite3.connect(DB_FILENAME)
    curs = conn.cursor()
    curs.execute('''UPDATE Protein SET GeneSymbol = ? WHERE UniProtAccNum = ?''',
        (gene_symbol, uniprot_accnum))
    conn.commit()
    conn.close()

def update_protein_process_function(uniprot_accnum, process_function):
    conn = sqlite3.connect(DB_FILENAME)
    curs = conn.cursor()
    curs.execute('''UPDATE Protein SET ProcessFunction = ? WHERE UniProtAccNum = ?''',
        (process_function, uniprot_accnum))
    conn.commit()
    conn.close()

def update_protein_ions(uniprot_accnum, ions):
    conn = sqlite3.connect(DB_FILENAME)
    curs = conn.cursor()
    curs.execute('''UPDATE Protein SET Ions = ? WHERE UniProtAccNum = ?''',
        (ions, uniprot_accnum))
    conn.commit()
    conn.close()

def update_protein_gating(uniprot_accnum, gating):
    conn = sqlite3.connect(DB_FILENAME)
    curs = conn.cursor()
    curs.execute('''UPDATE Protein SET Gating = ? WHERE UniProtAccNum = ?''',
        (gating, uniprot_accnum))
    conn.commit()
    conn.close()

def update_protein_in_betse(uniprot_accnum, in_betse):
    conn = sqlite3.connect(DB_FILENAME)
    curs = conn.cursor()
    curs.execute('''UPDATE Protein SET InBETSE = ? WHERE UniProtAccNum = ?''',
        (in_betse, uniprot_accnum))
    conn.commit()
    conn.close()

def update_protein_ion_channel_class_desc(uniprot_accnum, ion_channel_class_desc):
    conn = sqlite3.connect(DB_FILENAME)
    curs = conn.cursor()
    curs.execute('''UPDATE Protein SET IonChannelClassDesc = ? WHERE UniProtAccNum = ?''',
        (ion_channel_class_desc, uniprot_accnum))
    conn.commit()
    conn.close()

def update_protein_ion_channel_sub_class(uniprot_accnum, ion_channel_sub_class):
    conn = sqlite3.connect(DB_FILENAME)
    curs = conn.cursor()
    curs.execute('''UPDATE Protein SET IonChannelSubClass = ? WHERE UniProtAccNum = ?''',
        (ion_channel_sub_class, uniprot_accnum))
    conn.commit()
    conn.close()

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

def lookup_uniprot_accnum_by_gene_symbol_multiple(gene_symbol):
    conn = sqlite3.connect(DB_FILENAME)
    curs = conn.cursor()
    curs.execute('''SELECT UniProtAccNum FROM Protein WHERE GeneSymbol=?''', 
        (gene_symbol,))
    resultset = curs.fetchall()
    if resultset == None:
        result = []
    else:
        result = resultset
    conn.close()
    return resultset

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

def lookup_compound_by_chembl_id(chembl_id):
    conn = sqlite3.connect(DB_FILENAME)
    curs = conn.cursor()
    curs.execute('''SELECT Id FROM Compound WHERE ChemblId=?''', 
        (chembl_id,))
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


    curs.execute('''CREATE TABLE Compound(
        Id INTEGER PRIMARY KEY AUTOINCREMENT NOT NULL,
        SMILES TEXT NOT NULL,
        InChI TEXT NOT NULL,
        Name TEXT NOT NULL,
        ChemblId TEXT NOT NULL,
        Synonyms TEXT NOT NULL,
        ApprovalStatus TEXT NOT NULL,
        FirstApprovalYear TEXT NOT NULL,
        SourceDBName NOT NULL
        );''')


def lookup_compound_by_id(id):
    conn = sqlite3.connect(DB_FILENAME)
    curs = conn.cursor()
    curs.execute('''SELECT Id, SMILES, InChI, Name, ChemblId, Synonyms,
            ApprovalStatus, FirstApprovalYear, SourceDBName
            FROM Compound WHERE Id=?''', 
        (id,))
    resultset = curs.fetchone()
    if resultset == None:
        result = ''
    else:
        result = resultset
    conn.close()
    return result

def lookup_interaction_by_uniprot(upac):
    conn = sqlite3.connect(DB_FILENAME)
    curs = conn.cursor()
    curs.execute('''SELECT Id FROM Interaction WHERE TargetUniprotAccNum=?''', 
        (upac,))
    resultset = curs.fetchall()
    result = resultset
    conn.close()
    return result

def setup_externaldb():
    add_externaldb('biogps', 'http://biogps.org')
    add_externaldb('chembl', 'https://www.ebi.ac.uk/chembl/')
    add_externaldb('brenda', 'http://www.brenda-enzymes.org')
    add_externaldb('drugbank', 'https://www.drugbank.ca')
    add_externaldb('zinc15', 'http://zinc15.docking.org')
    add_externaldb('hmdb', 'http://www.hmdb.ca')
    add_externaldb('channelpedia', 'http://channelpedia.epfl.ch')

def cleanup_externaldb():
    conn = sqlite3.connect(DB_FILENAME)
    curs = conn.cursor()
    curs.execute('''DELETE FROM ExternalDB''')
    conn.commit()
    conn.close()

def setup_ontology():
    tissue_dict = {}
    tissue_file = open(TISSUE_FILENAME)
    for line in tissue_file:
        tissue_dict[line.strip()] = 1
    tissue_file.close()
    parser = go_obo_parser.parseGOOBO(OBO_FILENAME)
    for record in parser:
        name = record['name']
        if name in tissue_dict:
            bto_id = record['id']
            add_tissue(name, bto_id=bto_id)

def cleanup_ontology():
    conn = sqlite3.connect(DB_FILENAME)
    curs = conn.cursor()
    curs.execute('''DELETE FROM Tissue''')
    conn.commit()
    conn.close()

def setup_biogps():
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

def remove_nonascii(s):
    return s.encode('ascii', errors='ignore')

def setup_ion_channels():
    prot_file = open(PROTEIN_FILENAME, 'rU')
    #prot_reader = unicodecsv.reader(prot_file, encoding='utf-8')
    prot_reader = csv.reader(prot_file)
    prot_reader.next()
    for row in prot_reader:
        name = row[2]
        #process_function = remove_nonascii(row[3])
        process_function = ''
        #if process_function.startswith('NOT AN ION CHANNEL'):
        #    continue
        uniprot_accnum = row[0]
        if uniprot_accnum == '':
            continue
        #if not uniprot_accnum[0] in ['O', 'P', 'Q']:
        #    continue
        gene_symbol = row[1]
        if gene_symbol == '-':
            gene_symbol = ''
        #if gene_symbol == '':
        #    continue
        ions = row[4]
        gating = row[5]
        ion_channel_subclass = row[6]
        in_betse = row[7]
        #check_uniprot_accnum = lookup_uniprot_accnum_by_gene_symbol(gene_symbol)
        #if check_uniprot_accnum != '':
        #    print 'Dup:', gene_symbol, check_uniprot_accnum
        #check_gene_symbol = lookup_gene_symbol_by_uniprot_accnum(uniprot_accnum)
        #if check_gene_symbol != '':
        #    print 'Dup:', uniprot_accnum, check_gene_symbol
        
        if not exists_protein(uniprot_accnum):
            add_protein(uniprot_accnum, gene_symbol, name, process_function, ions,
                    gating, in_betse, ion_channel_subclass)
            pass
        
        else:
            old_protein_record = lookup_protein_by_accnum(uniprot_accnum)
            old_uniprot_accnum = old_protein_record[0]
            old_gene_symbol = old_protein_record[1]
            old_name = old_protein_record[2]
            old_process_function = old_protein_record[3]
            old_ions = old_protein_record[4]
            old_gating = old_protein_record[5]
            old_in_betse = old_protein_record[6]
            old_ion_channel_class_desc = old_protein_record[7]
            old_ion_channel_subclass = old_protein_record[8]
            old_chembl_id = old_protein_record[9]
            
            update_protein_gene_symbol(uniprot_accnum, gene_symbol)
            update_protein_name(uniprot_accnum, name)
            #update_protein_process_function(uniprot_accnum, process_function)
            update_protein_ions(uniprot_accnum, ions)
            update_protein_gating(uniprot_accnum, gating)
            if old_in_betse == '':
                if in_betse == 'yes':
                    in_betse = 'Y'
                elif in_betse == 'soon':
                    in_betse = '?'
                elif in_betse == 'doable':
                    in_betse = '?'
                else:
                    in_betse = 'N'
                update_protein_in_betse(uniprot_accnum, in_betse)
            update_protein_ion_channel_sub_class(uniprot_accnum, ion_channel_subclass)
        
        #print 'Added %s' % uniprot_accnum
    prot_file.close()

def setup_channel_classes():
    classfile = open(CHANNEL_CLASSES_FILENAME, 'rU')
    classreader = csv.reader(classfile)
    classreader.next()
    for row in classreader:
        channel_superclass = row[0]
        channel_class = row[1]
        channel_subclass = row[2]
        subfamily = row[3]
        gene = row[4]
        if not exists_channel_superclass(channel_superclass):
            add_channel_superclass(channel_superclass)
        if not exists_channel_class(channel_class):
            add_channel_class(channel_class, channel_superclass)
        if not exists_channel_subclass(channel_subclass):
            add_channel_subclass(channel_subclass, channel_class)
        update_protein_subclass(gene, channel_subclass)
    classfile.close()

def setup_all_proteins():
    mapping_dict = {}
    uniprot_symbol_file = open(HUMAN_IDMAPPING)
    for line in uniprot_symbol_file:
        uniprot_accum = line.strip().split()[0]
        gene_symbol = line.strip().split()[2]
        mapping_dict[gene_symbol] = uniprot_accum
    uniprot_symbol_file.close()
    prot_file = open(BIOGPS_GNF1H_ANNOT, 'rU')
    prot_reader = csv.reader(prot_file, delimiter='\t')
    prot_reader.next()
    for row in prot_reader:
        gene_symbol = row[6]
        if not gene_symbol in mapping_dict:
            'Could not find uniprot accession for %s' % gene_symbol
            continue
        uniprot_accum = mapping_dict[gene_symbol]
        if exists_protein(uniprot_accum):
            continue
        add_protein(uniprot_accum, gene_symbol, '', '')
        #print 'Added %s %s' % (uniprot_accum, gene_symbol)
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
        chembl_compound_id = line[114:132].strip()
        compound_name_in_doc = line[132:334].strip()
        first_approval = line[334:342].strip()
        withdrawn_flag = line[342:346].strip()
        inchi = line[346:].strip()
        # This is a temporary way of checking if the compound exists
        if exists_protein(uniprot) and not exists_compound(chembl_compound_id):
            if withdrawn_flag == 1:
                approval_status = 'WITHDRAWN'
            else:
                approval_status = 'UNKNOWN'
            add_compound('', '', compound_name_in_doc, chembl_compound_id, '',
                    approval_status, first_approval, 'chembl')
            #print 'Added compound %s' % chembl_compound_id
        if exists_protein(uniprot) and exists_compound(chembl_compound_id):
            compound_id = lookup_compound_by_name(chembl_compound_id)
            print compound_id
            if not exists_interaction(uniprot, compound_id):
                add_interaction(uniprot, compound_id, sourcedb_name='chembl')
                #print 'Added interaction %s %s' % (uniprot, chembl_compound_id)
    chembl_compound_file.close()

def setup_hmdb():
    pass

def setup_drugbank():
    pass

def setup_target_compound():
    tcf = open(TARGET_COMPOUND_FILENAME)
    tcr = csv.reader(tcf)
    tcr.next()
    for row in tcr:
        uniprot = row[1]
        effect_type = row[7]
        compound_name = row[9]
        if compound_name == '':
            continue
        chembl_id = row[10]
        param = row[11]
        value = float(row[12])
        if exists_protein(uniprot):
            if not exists_compound(chembl_id):
                add_compound('', '', compound_name, chembl_id, 
                        '', '', '', 'chembl')
            compound_id = lookup_compound_by_chembl_id(chembl_id)
            add_interaction(uniprot, compound_id, action_type=effect_type,
                    strength=value, strength_units=param, sourcedb_name='chembl')
    tcf.close()

def cleanup_compounds():
    conn = sqlite3.connect(DB_FILENAME)
    curs = conn.cursor()
    curs.execute('''DELETE FROM Compound''')
    conn.commit()
    conn.close()

def cleanup_interactions():
    conn = sqlite3.connect(DB_FILENAME)
    curs = conn.cursor()
    curs.execute('''DELETE FROM Interaction''')
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
        uniprot_accnum_list = lookup_uniprot_accnum_by_gene_symbol_multiple(symbol.strip())
        if len(uniprot_accnum_list) == 0:
            continue
        for i in xrange(0, len(tissue_list)):
            expr_level = float(row[i + 1])
            tissue_name = lookup_dbtissue('biogps', tissue_list[i])
            for uniprot_accnum in uniprot_accnum_list:
                #if uniprot_accnum in ['Q12809', 'Q6U279', 'Q6U283', 'Q6U287', 'Q86U57']:
                #    print uniprot_accnum, expr_level
                uniprot_accnum = uniprot_accnum[0]
                add_expression(tissue_name, uniprot_accnum, expr_level, 
                        assay_type='microarray', dataset_name=BIOGPS_DATASET, 
                        sourcedb_name='biogps')
                #print 'Added %s %s' % (uniprot_accnum, tissue_name)
    microarray_f.close()

def process_biogps_data():
    annot_dict = process_biogps_annotation()
    process_biogps_microarray(annot_dict)


def process_protein_atlas_data():
    ensembl_to_uniprot = {}
    fmap = open(UNIPROT_ENSEMBL_MAPPING)
    for line in fmap:
        row = line.strip().split('\t')
        upac = row[0]
        ensg = row[2]
        ensembl_to_uniprot[ensg] = upac
    fmap.close()
    
    tissue_to_bto = {}
    ftissue = open(HPA_TISSUE_TRANSLATION_TABLE)
    rtissue = csv.reader(ftissue)
    rtissue.next()
    for row in rtissue:
        hpa_tissue = row[0]
        bto_tissue = row[1]
        tissue_to_bto[hpa_tissue] = bto_tissue
    ftissue.close()
    
    fhpa = open(HPA_RNA_TISSUE)
    fhpa.next()
    for line in fhpa:
        row = line.strip().split('\t')
        if len(row) < 5 or row[0].strip() == '':
            continue
        ensg = row[0]
        if not ensg in ensembl_to_uniprot:
            continue
        upac = ensembl_to_uniprot[ensg]
        gene_name = row[1]
        tissue = row[2]
        tissue_bto = tissue_to_bto[tissue]
        value = row[3]
        if exists_protein(upac):
            add_expression(tissue_bto, upac, value, assay_type='RNA-seq',
                    dataset_name=HPA_DATASET, sourcedb_name=EXTDB_HPA)
    fhpa.close()

def cleanup_expression():
    conn = sqlite3.connect(DB_FILENAME)
    curs = conn.cursor()
    curs.execute('''DELETE FROM Expression''')
    conn.commit()
    conn.close()

def get_expression_for_tissue(tissue_name):
    prot_expr_list = []
    conn = sqlite3.connect(DB_FILENAME)
    curs = conn.cursor()
    for row in curs.execute('''SELECT ProteinUniProtAccNum, ExprLevel
            FROM Expression WHERE TissueName = ?''',
            (tissue_name,)):
        prot_expr_list.append((row[0], row[1]))
    conn.commit()
    conn.close()
    return prot_expr_list

def normalize_prot_expr_list(prot_expr_list):
    normalized_list = []
    total_expr = 0.0
    for elem in prot_expr_list:
        upac = elem[0]
        expr = elem[1]
        total_expr += expr
    for elem in prot_expr_list:
        upac = elem[0]
        expr = elem[1]
        norm_expr = expr / total_expr
        normalized_list.append((upac, norm_expr))
    return normalized_list

def get_uniprot_seq(upac):
    ls = os.listdir('data/seq')
    fn = '%s.fasta' % upac
    if fn == 'Q96LR1.fasta':
        return ''
    if not fn in ls:
        cmd = 'wget -O data/seq/%s http://www.uniprot.org/uniprot/%s' % (fn, fn)
        #print cmd
        os.system(cmd)
    f = open('data/seq/%s' % fn)
    f.next()
    seqlist = []
    for line in f:
        seqlist.append(line.strip())
    f.close()
    seq = ''.join(seqlist)
    return seq

def get_expression_prob_for_all_tissues():
    aa_list = ['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P',
            'Q', 'R', 'S', 'T', 'V', 'W', 'Y']
    out_file = open('aa_out.csv', 'w')
    out_file.write('tissue')
    for aa in aa_list:
        out_file.write(',%s' % aa)
    out_file.write('\n')
    tissue_list = []
    conn = sqlite3.connect(DB_FILENAME)
    curs = conn.cursor()
    for row in curs.execute('''SELECT Name FROM Tissue'''):
        tissue_list.append(row[0])
    conn.commit()
    conn.close()
    data_dict = {}
    count = 0
    for tissue in tissue_list:
        #if count == 1:
        #    break
        prot_expr_list = get_expression_for_tissue(tissue)
        if len(prot_expr_list) > 0:
            print 'Normalizing %s' % tissue
            data_dict[tissue] = normalize_prot_expr_list(prot_expr_list)
            count += 1
    for tissue in data_dict:
        print 'Probablizing %s' % tissue
        overall_aa_prob = {}
        normalized_list = data_dict[tissue]
        for elem in normalized_list:
            upac, expr = elem
            seq = get_uniprot_seq(upac)
            seq_freq_dict = {}
            total_aa = 0.0
            for aa in seq:
                total_aa += 1
                if aa in seq_freq_dict:
                    seq_freq_dict[aa] += 1.0
                else:
                    seq_freq_dict[aa] = 1.0
            seq_prob_dict = {}
            for aa in seq_freq_dict:
                prob = seq_freq_dict[aa] / total_aa
                seq_prob_dict[aa] = prob
                if aa in overall_aa_prob:
                    overall_aa_prob[aa] += expr * prob
                else:
                    overall_aa_prob[aa] = expr * prob
        out_file.write('%s' % tissue)
        for aa in aa_list:
            out_file.write(',%.6f' % overall_aa_prob[aa])
        out_file.write('\n')
    out_file.close()

def write_db_stats():

    conn = sqlite3.connect(DB_FILENAME)
    curs = conn.cursor()
    curs.execute("SELECT * FROM sqlite_master WHERE type='table'")
    table_name_list = []
    for table_info in curs:
        table_name = table_info[1]
        table_name_list.append(table_name)
    conn.close()
    
    for table_name in table_name_list:
        conn = sqlite3.connect(DB_FILENAME)
        curs = conn.cursor()
        curs.execute('SELECT COUNT(*) FROM %s' % (table_name,))
        row_count = curs.fetchone()[0]
        print table_name, row_count
        conn.close()
        
    useful_data_points = 0
    uniprot_list = get_all_protein_uniprot()
    for uniprot in uniprot_list:
        uniprot = uniprot[0]
        interaction_id_list = lookup_interaction_ids_by_uniprot(uniprot)
        expression_id_list = lookup_expression_ids_by_uniprot(uniprot)
        if len(interaction_id_list) > 0 and len(expression_id_list) > 0:
            useful_data_points += 1
    print 'Protein with interaction and expression %d' % useful_data_points

    #uniprot_list = get_all_protein_uniprot()
    #for uniprot in uniprot_list:
    #    uniprot = uniprot[0]
    #    gene_symbol = lookup_gene_symbol_by_uniprot_accnum(uniprot)
    #    print uniprot, gene_symbol


def input_drugbank():
    drug_info_file = open(OUTPUT_ORGANIZED_DRUG_INFO_DRUGBANK)
    drug_info_file.next()
    for line in drug_info_file:
        drugbank_id = line[0:13].strip()
        drug_name = line[13:64].strip()
        status = line[64:121].strip()
        status_list = status.split(',')
        molecular_formula = line[121:152].strip()
        if molecular_formula == 'NONE':
            molecular_formula = ''
        inchi_string = line[152:].strip()
        if inchi_string == 'NONE':
            inchi_string = ''
        record = (drugbank_id, drug_name, status_list, molecular_formula,
                inchi_string)
    drug_info_file.close()
    
    target_info_file = open(OUTPUT_ORGANIZED_TARGET_INFO_DRUGBANK)
    target_info_file.next()
    for line in target_info_file:
        drugbank_id = line[0:13].strip()
        target_id = line[13:26].strip()
        category = line[26:39].strip()
        effect = line[39:80].strip()
        if effect == 'UNKNOWN':
            effect = ''
        uniprot_id = line[80:93].strip()
        if uniprot_id == 'NONE':
            uniprot_id = ''
        target_name = line[93:].strip()
        record = (drugbank_id, target_id, category, effect, uniprot_id, 
                target_name)
    target_info_file.close()

def input_hmdb():
    drug_info_file = open(OUTPUT_ORGANIZED_DRUG_INFO_HMDB)
    drug_info_file.next()
    for line in drug_info_file:
        hmdb_id = line[0:13].strip()
        name = line[13:124].strip()
        molecular_formula = line[124:155].strip()
        inchi_string = line[155:].strip()
        record = (hmdb_id, name, molecular_formula, inchi_string)
    drug_info_file.close()

    target_info_file = open(OUTPUT_ORGANIZED_TARGET_INFO_HMDB)
    target_info_file.next()
    for line in target_info_file:
        hmdb_id = line[0:13].strip()
        target_id = line[13:26].strip()
        uniprot_id = line[26:39].strip()
        target_name = line[39:].strip()
        record = (hmdb_id, target_id, uniprot_id, target_name)
    target_info_file.close()

def input_ttd():
    ttd_file = open(OUTPUT_ORGANIZED_TTD)
    ttd_file.next()
    for line in ttd_file:
        target_id = line[0:13].strip()
        uniprot_id = line[13:26].strip()
        if uniprot_id == 'UniProt-NONE':
            uniprot_id = ''
        target_name = line[26:82].strip()
        drug_type = line[82:95].strip()
        chembl_id = line[95:111].strip()
        if chembl_id == 'UNKNOWN':
            chembl_id = ''
        ttd_drug_id = line[111:124].strip()
        drug_name = line[124:210].strip()
        inchi = line[210:].strip()
        if inchi == 'InChI-NONE':
            inchi = ''
        record = (target_id, uniprot_id, target_name, drug_type, chembl_id,
                ttd_drug_id, drug_name, inchi)
    ttd_file.close()

def input_chembl_uniprot():
    compound_file = open(OUTPUT_QUERY_CHEMBL_UNIPROT_COMPOUND)
    for line in compound_file:
        accession = line[0:14].strip()
        pref_name = line[14:96].strip()
        target_chembl_id = line[96:114].strip()
        
        # Is this a new target?
        # If yes, add it to the database
        # If no, check values match, updating if necessary
        
        if exists_protein(accession):
            protein_record = lookup_protein_by_accnum(accession)
            if protein_record[2] == '':
                update_protein_name(accession, pref_name)
            if protein_record[9] == '':
                update_protein_chembl_id(accession, target_chembl_id)
        else:
            add_protein(accession, name=pref_name)
        
        compound_chembl_id = line[114:132].strip()
        compound_name = line[132:334].strip()
        first_approval = line[334:342].strip()
        withdrawn_flag = line[342:346].strip()
        standard_inchi = line[346:].strip()
        
        # Is this a new compound?
        # If yes, add it to the database
        # If no, check values match, updating if necessary
        
        if exists_compound(compound_chembl_id):
            compound_id = lookup_compound_by_chembl_id(compound_chembl_id)
            compound_record = lookup_compound_by_id(compound_id)
        else:
            if withdrawn_flag == '1':
                approval_status = 'WITHDRAWN'
            else:
                approval_status = ''
            add_compound('', standard_inchi, compound_name,
                    chembl_id=compound_chembl_id,
                    first_approval_year=first_approval,
                    approval_status=approval_status, sourcedb_name=EXTDB_CHEMBL)

        # Add an interaction record for the (target, compound) pair
        
        # INCOMPLETE
        
    compound_file.close()
    drug_file = open(OUTPUT_QUERY_CHEMBL_UNIPROT_DRUG)
    for line in drug_file:
        accession = line[0:14].strip()
        pref_name = line[14:66].strip()
        target_chembl_id = line[66:84].strip()
        
        # Is this a new target?
        # If yes, add it to the database
        # If no, check values match, updating if necessary
        
        # INCOMPLETE
        
        compound_chembl_id = line[114:132].strip()
        compound_chembl_id = line[84:102].strip()
        compound_name = line[102:134].strip()
        first_approval = line[134:141].strip()
        withdrawn_flag = line[141:145].strip()
        mechanism_of_action = line[145:197].strip()
        action_type = line[197:229].strip()
        standard_inchi = line[229:].strip()
        
        # Is this a new compound?
        # If yes, add it to the database
        # If no, check values match, updating if necessary
        
        # INCOMPLETE
        
        # Add an interaction record for the (target, compound) pair
        
        # INCOMPLETE
        
    drug_file.close()

def input_chembl_assays():
    compound_file = open(OUTPUT_COMPOUND_ASSAYS)
    for line in compound_file:
        line_split = line.split('\t')
        uniprot = line_split[0].strip()
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
        assay_type = line_split[12].strip() # Not used
        assay_description = line_split[13].strip() # Not used
        
        if not exists_compound(compound_chembl_id):
            #add_compound('', '', compound_name.decode('utf_8').encode('utf_8'), 
            #        chembl_id=compound_chembl_id.decode('utf_8').encode('utf_8'),
            #        synonyms=iupac_name.decode('utf_8').encode('utf_8'),
            #        sourcedb_name=EXTDB_CHEMBL)
            add_compound('', '', compound_name, 
                    chembl_id=compound_chembl_id,
                    synonyms=iupac_name,
                    sourcedb_name=EXTDB_CHEMBL)
        
        if not exists_protein(uniprot):
            add_protein(uniprot, name=target_name,
                    process_function=target_type, chembl_id=target_chembl_id)
    
        compound_id = lookup_compound_by_chembl_id(compound_chembl_id)

        if not exists_interaction(uniprot, compound_id):
            #if uniprot in ['Q12809', 'Q6U279', 'Q6U283', 'Q6U287', 'Q86U57']:
            #    print uniprot, assay_value

            #add_interaction(uniprot, compound_id,
            #        action_type=target_type.decode('utf_8').encode('utf_8'),
            #        strength=assay_value.decode('utf_8').encode('utf_8'),
            #        strength_units=assay_units.decode('utf_8').encode('utf_8'),
            #        assay_type=assay_standard_type.decode('utf_8').encode('utf_8'),
            #        chembl_id=assay_chembl_id,
            #        sourcedb_name=EXTDB_CHEMBL)

            add_interaction(uniprot, compound_id,
                    action_type='',
                    strength=assay_value.decode('utf_8').encode('utf_8'),
                    strength_units=assay_units.decode('utf_8').encode('utf_8'),
                    assay_type=assay_standard_type.decode('utf_8').encode('utf_8'),
                    chembl_id=assay_chembl_id,
                    sourcedb_name=EXTDB_CHEMBL)

                    
        if assay_standard_relation != '=':
            print 'Assay standard relation is not ='

    compound_file.close()
    drug_file = open(OUTPUT_DRUG_ASSAYS)
    for line in drug_file:
        line_split = line.split('\t')
        uniprot = line_split[0].strip()
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

        if not exists_compound(active_compound_chembl_id):
            add_compound('', '', active_compound_name, 
                    chembl_id=active_compound_chembl_id,
                    synonyms=dosed_compound_name,
                    sourcedb_name=EXTDB_CHEMBL)
        
        if not exists_protein(uniprot):
            add_protein(uniprot, name=target_name,
                    process_function=target_type, chembl_id=target_chembl_id)

        compound_id = lookup_compound_by_chembl_id(active_compound_chembl_id)

        if not exists_interaction(uniprot, compound_id):
            add_interaction(uniprot, compound_id,
                    action_type=drug_mechanism.decode('utf_8').encode('utf_8'),
                    strength=assay_value.decode('utf_8').encode('utf_8'),
                    strength_units=assay_units.decode('utf_8').encode('utf_8'),
                    assay_type=assay_standard_type.decode('utf_8').encode('utf_8'),
                    chembl_id=assay_chembl_id,
                    sourcedb_name=EXTDB_CHEMBL)

    drug_file.close()

def write_important_data():

    tissue_name_list = []
    tissue_record_list = get_all_tissues()
    for tissue in tissue_record_list:
        tissue_name = tissue[0]
        tissue_name_list.append(tissue_name)
    protein_upac_list = []
    protein_record_list = get_all_protein_uniprot()
    for protein_record in protein_record_list:
        upac = protein_record[0]
        has_expression = exists_expression(upac)
        has_interaction = exists_any_interaction(upac)
        if has_expression and has_interaction:
            protein_upac_list.append(upac)
    line = 'Uniprot Accession,Gene Symbol,Number of interacting compounds,In BETSE?,'
    for tissue_name in tissue_name_list:
        line += tissue_name
        line += ','
    print line
    for upac in protein_upac_list:
        line = ''
        line += upac
        line += ','
        line += lookup_gene_symbol_by_uniprot_accnum(upac)
        line += ','
        number_of_interacting_compounds = len(lookup_interaction_by_uniprot(upac))
        line += '%d' % number_of_interacting_compounds
        line += ','
        in_betse_yn = is_in_betse(upac)
        line += in_betse_yn
        line += ','
        for tissue_name in tissue_name_list:
            expr_level_list = lookup_expression(upac, tissue_name)
            total_expr_level = 0.0
            for expr_level in expr_level_list:
                total_expr_level += float(expr_level[0])
            avg_expr_level = total_expr_level / float(len(expr_level_list))
            line += '%.2f' % avg_expr_level
            line += ','
        print line
        


#expressed_protein_dict = {}
#tissue_list = get_all_tissues()
#for tissue in tissue_list:
#    tissue_name = tissue[0]
#    expression_record_list = get_all_expression()
#    for expression_record in expression_record_list:
#        expression_id = expression_record[0]
#        expression_tissue_name = expression_record[1]
#        expression_protein_uniprot_accnum = expression_record[2]
#        if expression_tissue_name == tissue_name:
#            if not expression_protein_uniprot_accnum in expressed_protein_dict:
#                if exists_any_interaction(expression_protein_uniprot_accnum):
#                    expressed_protein_dict[expression_protein_uniprot_accnum] = 1
#for protein in expressed_protein_list:
#    line = ''
#    line += expression_protein_uniprot_accnum + ','
#    expression_record_list = get_all_expression()
#    for expression_record in expression_record_list:
#        expression_id = expression_record[0]
#        expression_tissue_name = expression_record[1]
#        expression_protein_uniprot_accnum = expression_record[2]
#        expression_level = expression_record[3]
#        line += expression_level + ','
#    print line

#   myresults = form.cleaned_data
#    tissue_selection = myresults['tissue_mc']
#    tissue_selection.sort()
#    protein_selection = myresults['protein_mc']
#
#    protein_label_list = []
#    protein_label_dict = {}
#    for upac in protein_selection:
#        protein = Protein.objects.get(uniprotaccnum=upac)
#        ic_subclass = protein.ionchannelsubclass
#        #ic_class = Channelsubclass.objects.get(name=ic_subclass).class_field
#        #ic_superclass = Channelclass.objects.get(name=ic_class).superclass
#        protein_label = '%s (%s)' % (protein.uniprotaccnum, protein.genesymbol)
#        protein_label_list.append(protein_label)
#        protein_label_dict[protein_label] = upac
#    protein_label_list.sort()
#
#    final_results = []
#    for tissue_name in tissue_selection:
#        for protein_label in protein_label_list:
#            upac = protein_label_dict[protein_label]
#            expression_level_qs = Expression.objects.filter(tissuename=tissue_name).filter(proteinuniprotaccnum=upac)
#            if len(expression_level_qs) < 1:
#                expression_level = '0.00'
#            else:
#                expression_level = expression_level_qs[0].exprlevel
#            interaction_qs = Interaction.objects.filter(targetuniprotaccnum=upac)
#            if len(interaction_qs) > 0:
#                for interaction in interaction_qs:
#                    compound_id = interaction.compoundid
#                    compound_record = Compound.objects.get(id=compound_id)
#                    compound_name = compound_record.name
#                    compound_chemblid = compound_record.chemblid
#                    compound_label = '%s (%s)' % (compound_name, compound_chemblid)
#                    effect_type = ''
#                    param_type = interaction.assaytype
#                    param_value = '%.2f' % interaction.strength
#
#                    final_results.append([tissue_name, protein_label, expression_level,
#                            compound_label, effect_type, param_type, param_value])
#            else:
#                final_results.append([tissue_name, protein_label, expression_level,
#                            '', '', '', ''])




