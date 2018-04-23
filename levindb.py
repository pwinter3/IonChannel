import codecs
import csv
import go_obo_parser
import os
import sqlite3
import unicodecsv
import time
import sys

# Constants to be stored in db

DATASET_BIOGPS = 'U133AGNF1B'
DATASET_HPA = 'hpa'
ASSAY_MICROARRAY = 'microarray'
ASSAY_RNASEQ = 'RNA-seq'
ASSAY_IMMUNO = 'immuno'
EXTDB_BIOGPS = 'biogps'
EXTDB_BRENDA = 'brenda'
EXTDB_CHANNELPEDIA = 'channelpedia'
EXTDB_CHEMBL = 'chembl'
EXTDB_DRUGBANK = 'drugbank'
EXTDB_HMDB = 'hmdb'
EXTDB_HPA = 'hpa'
EXTDB_ZINC15 = 'zinc15'

# Constant paths

PATH_DB = 'levin.db'
PATH_BIOGPS_GCRMA = 'data/biogps/U133AGNF1B.gcrma.avg.csv'
PATH_BIOGPS_GNF1H_ANNOT = 'data/biogps/gnf1h.annot2007.tsv'
PATH_BIOGPS_TRANSLATION_TABLE = 'data/biogps/biogps_translation_table.csv'
PATH_BTO_BRENDA_TISSUE_OBO = 'data/bto/BrendaTissueOBO.obo'
PATH_BTO_TISSUES = 'data/bto/tissues.txt'
PATH_CHANNEL_CLASSES = 'data/channel-classes/channel-classes.csv'
PATH_CHEMBL_ASSAYS_COMPOUND = 'data/chembl-assays/output_compound-assays_v1.1.dat'
PATH_CHEMBL_ASSAYS_DRUG = 'data/chembl-assays/output_drug-assays_v1.1.dat'
PATH_CHEMBL_HUMAN_ASSAYS_COMPOUND = 'data/chembl-assays/output_human_compound-assays_v1.1.dat'
PATH_CHEMBL_HUMAN_ASSAYS_DRUG = 'data/chembl-assays/output_human_drug-assays_v1.1.dat'
PATH_GO_ION_CHANNELS = 'data/go/ion-channels-1.csv'
PATH_GO_QUICKGO = 'data/go/QuickGO-ion-channel-COMBINED-human.dat'
PATH_HPA_RNA_TISSUE = 'data/hpa/rna_tissue.tsv'
PATH_HPA_TRANSLATION_TABLE = 'data/hpa/hpa_tissue_translation_table.csv'
PATH_TARGET_COMPOUND = 'data/target-compound/target-compound.csv'
PATH_TC_CHEMBL_COMPOUND = 'data/target-compound/output-query_ChEMBL-uniprot-compound.dat'
PATH_TC_CHEMBL_DRUG = 'data/target-compound/output-query_ChEMBL-uniprot-drug.dat'
PATH_TC_DRUG_INFO_DRUGBANK = 'data/target-compound/output-organized_drug-info-Drugbank.dat'
PATH_TC_DRUG_INFO_HMDB = 'data/target-compound/output-organized_drug-info-HMDB.dat'
PATH_TC_TARGET_INFO_DRUGBANK = 'data/target-compound/output-organized_target-info-Drugbank.dat'
PATH_TC_TARGET_INFO_HMDB = 'data/target-compound/output-organized_target-info-HMDB.dat'
PATH_TC_TTD = 'data/target-compound/output-organized-TTD.dat'
PATH_UNIPROT_ENSEMBL_IDMAPPING = 'data/uniprot/HUMAN_9606_idmapping_Ensembl.dat'
PATH_UNIPROT_GENBANK = 'data/uniprot/GenBankUniProt.txt' 
PATH_UNIPROT_HUMAN_IDMAPPING = 'data/uniprot/HUMAN_9606_idmapping.dat'
PATH_HPA_PATHOLOGY = 'data/hpa/pathology.tsv'
PATH_CHANNELPEDIA_INFO = 'data/channel-classes/channelpedia_info.csv'
PATH_HPA_NORMAL = 'data/hpa/normal_tissue.tsv'

# Table definitions

def create_tables():
    conn = sqlite3.connect(PATH_DB)
    curs = conn.cursor()
    curs.execute('''
        CREATE TABLE ChannelSuperClass(
        Name TEXT PRIMARY KEY NOT NULL);
        ''')
    curs.execute('''
        CREATE TABLE ChannelClass(
        Name TEXT PRIMARY KEY NOT NULL,
        SuperClass TEXT NOT NULL);
        ''')
    curs.execute('''
        CREATE TABLE ChannelSubClass(
        Name TEXT PRIMARY KEY NOT NULL,
        Class TEXT NOT NULL,
        ChannelpediaText TEXT NOT NULL,
        ChannelpediaURL TEXT NOT NULL);
        ''')
    curs.execute('''
        CREATE TABLE Tissue(
        Name TEXT PRIMARY KEY NOT NULL,
        BTOId TEXT UNIQUE NOT NULL,
        ParentTissueName TEXT NOT NULL);
        ''')
    curs.execute('''
        CREATE TABLE Protein(
        UniProtAccNum TEXT PRIMARY KEY NOT NULL,
        GeneSymbol TEXT NOT NULL,
        Name TEXT NOT NULL,
        ProcessFunction TEXT NOT NULL,
        Ions TEXT NOT NULL,
        Gating TEXT NOT NULL,
        InBETSE TEXT NOT NULL,
        IonChannelClassDesc TEXT NOT NULL,
        IonChannelSubClass TEXT NOT NULL,
        ChemblId TEXT NOT NULL);
        ''')
    curs.execute('''
        CREATE TABLE PDB(
        Id INTEGER PRIMARY KEY AUTOINCREMENT NOT NULL,
        UniProtAccNum TEXT NOT NULL,
        PDBID TEXT NOT NULL,
        HumanOnly TEXT NOT NULL);
        ''')
    curs.execute('''
        CREATE TABLE ExternalDB(
        Name TEXT PRIMARY KEY NOT NULL,
        URL TEXT NOT NULL);
        ''')
    curs.execute('''
        CREATE TABLE DBTissue(
        Id INTEGER PRIMARY KEY AUTOINCREMENT NOT NULL,
        ExternalDBName TEXT NOT NULL,
        TissueName TEXT NOT NULL,
        DBEquivalentTissueName TEXT NOT NULL);
        ''')
    curs.execute('''
        CREATE TABLE Expression(
        Id INTEGER PRIMARY KEY AUTOINCREMENT NOT NULL,
        TissueName TEXT NOT NULL,
        ProteinUniProtAccNum TEXT NOT NULL,
        ExprLevel REAL NOT NULL,
        ExprLevelQual TEXT NOT NULL,
        ExprUnits TEXT NOT NULL,
        AssayType TEXT NOT NULL,
        DatasetName TEXT NOT NULL,
        SourceDBName TEXT NOT NULL);
        ''')
    curs.execute('''
        CREATE TABLE DBGene(
        Id INTEGER PRIMARY KEY AUTOINCREMENT NOT NULL,
        ExternalDBName TEXT NOT NULL,
        GenbankAccNum TEXT NOT NULL,
        ProbeID TEXT NOT NULL);
        ''')
    curs.execute('''
        CREATE TABLE GenbankUniprot(
        Id INTEGER PRIMARY KEY AUTOINCREMENT NOT NULL,
        GenbankAccNum TEXT NOT NULL,
        UniProtAccNum TEXT NOT NULL);
        ''')
    curs.execute('''
        CREATE TABLE Compound(
        Id INTEGER PRIMARY KEY AUTOINCREMENT NOT NULL,
        SMILES TEXT NOT NULL,
        InChI TEXT NOT NULL,
        Name TEXT NOT NULL,
        ChemblId TEXT NOT NULL UNIQUE,
        Synonyms TEXT NOT NULL,
        ApprovalStatus TEXT NOT NULL,
        FirstApprovalYear TEXT NOT NULL,
        SourceDBName NOT NULL);
        ''')
    curs.execute('''
        CREATE TABLE Interaction(
        Id INTEGER PRIMARY KEY AUTOINCREMENT NOT NULL,
        TargetUniProtAccNum TEXT NOT NULL,
        CompoundId INTEGER NOT NULL,
        ActionType TEXT NOT NULL,
        ActionDesc TEXT NOT NULL,
        Strength REAL NOT NULL,
        StrengthUnits TEXT NOT NULL,
        AssayType TEXT NOT NULL,
        ChemblId TEXT NOT NULL,
        SourceDBName TEXT NOT NULL);
        ''')
    curs.execute('''
        CREATE TABLE Specificity(
        Id INTEGER PRIMARY KEY AUTOINCREMENT NOT NULL,
        TissueName TEXT NOT NULL,
        UniProtAccNum TEXT NOT NULL,
        SpecificityScore TEXT NOT NULL);
        ''')
    conn.commit()
    conn.close()

def vacuum():
    conn = sqlite3.connect(PATH_DB)
    curs = conn.cursor()
    curs.execute('''
        VACUUM;
        ''')
    conn.commit()
    conn.close()

# Routines for ChannelSuperClass table

def add_channel_super_class(name):
    conn = sqlite3.connect(PATH_DB)
    curs = conn.cursor()
    curs.execute('''
        INSERT INTO ChannelSuperClass
        (Name)
        VALUES (?)
        ''', (name,))
    conn.commit()
    conn.close()

def exists_channel_super_class(name):
    conn = sqlite3.connect(PATH_DB)
    curs = conn.cursor()
    curs.execute('''
        SELECT EXISTS(
            SELECT 1
            FROM ChannelSuperClass
            WHERE Name=? LIMIT 1
        )
        ''', (name,))
    row = curs.fetchone()
    if row[0] == 1:
        result = True
    else:
        result = False
    conn.close()
    return result

# Routines for ChannelClass table

def add_channel_class(name, super_class=''):
    conn = sqlite3.connect(PATH_DB)
    curs = conn.cursor()
    curs.execute('''
        INSERT INTO ChannelClass
        (Name, SuperClass)
        VALUES (?, ?)
        ''', (name, super_class))
    conn.commit()
    conn.close()

def exists_channel_class(name):
    conn = sqlite3.connect(PATH_DB)
    curs = conn.cursor()
    curs.execute('''
        SELECT EXISTS(
            SELECT 1
            FROM ChannelClass
            WHERE Name=?
            LIMIT 1
        )
        ''', (name,))
    row = curs.fetchone()
    if row[0] == 1:
        result = True
    else:
        result = False
    conn.close()
    return result

# Routines for ChannelSubClass table

def add_channel_sub_class(name, channel_class='', channelpedia_text='',
        channelpedia_url=''):
    conn = sqlite3.connect(PATH_DB)
    conn.text_factory = str
    curs = conn.cursor()
    curs.execute('''
        INSERT INTO ChannelSubClass
        (Name, Class, ChannelpediaText, ChannelpediaURL)
        VALUES (?, ?, ?, ?)
        ''', (name, channel_class, channelpedia_text, channelpedia_url))
    conn.commit()
    conn.close()

def exists_channel_sub_class(name):
    conn = sqlite3.connect(PATH_DB)
    curs = conn.cursor()
    expr_list = []
    curs.execute('''
        SELECT EXISTS(
            SELECT 1
            FROM ChannelSubClass
            WHERE Name=?
            LIMIT 1
        )
        ''', (name,))
    row = curs.fetchone()
    if row[0] == 1:
        result = True
    else:
        result = False
    conn.close()
    return result

# Routines for ExternalDB table

def add_externaldb(name, url=''):
    conn = sqlite3.connect(PATH_DB)
    curs = conn.cursor()
    curs.execute('''
        INSERT INTO ExternalDB
        (Name, URL)
        VALUES (?, ?)
        ''', (name, url))
    conn.commit()
    conn.close()

# Routines for Tissue table

def add_tissue(name, bto_id='', parent_tissue_name=''):
    conn = sqlite3.connect(PATH_DB)
    curs = conn.cursor()
    curs.execute('''
        INSERT INTO Tissue
        (Name, BTOId, ParentTissueName)
        VALUES (?, ?, ?)
        ''', (name, bto_id, parent_tissue_name))
    conn.commit()
    conn.close()

def get_all_tissue_names():
    conn = sqlite3.connect(PATH_DB)
    curs = conn.cursor()
    curs.execute('''
        SELECT Name
        FROM Tissue
        ''')
    resultset = curs.fetchall()
    conn.close()
    return resultset

def dump_tissues():
    conn = sqlite3.connect(PATH_DB)
    curs = conn.cursor()
    with codecs.open('tissue_dump.sql', 'w', 'utf-8') as f:
        for line in conn.iterdump():
            f.write('%s\n' % line)
    conn.close()

# Routines for Protein table

def add_protein(upac, gene_symbol='', name='', process_function='', ions='',
        gating='', in_betse='', ion_channel_class_desc='',
        ion_channel_sub_class='', chembl_id=''):
    conn = sqlite3.connect(PATH_DB)
    curs = conn.cursor()
    curs.execute('''
        INSERT INTO Protein
        (UniProtAccNum, GeneSymbol, Name, ProcessFunction, Ions, Gating,
        InBETSE, IonChannelClassDesc, IonChannelSubClass, ChemblId)
        VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
        ''', (upac, gene_symbol, name, process_function, ions, gating,
        in_betse, ion_channel_class_desc, ion_channel_sub_class, chembl_id))
    conn.commit()
    conn.close()

def update_protein_sub_class(upac, sub_class):
    conn = sqlite3.connect(PATH_DB)
    curs = conn.cursor()
    curs.execute('''
        UPDATE Protein
        SET IonChannelSubClass = ?
        WHERE UniProtAccNum = ?
        ''', (sub_class, upac))
    conn.commit()
    conn.close()

def update_protein_name(upac, name):
    conn = sqlite3.connect(PATH_DB)
    curs = conn.cursor()
    curs.execute('''
        UPDATE Protein
        SET Name = ?
        WHERE UniProtAccNum = ?
        ''', (name, upac))
    conn.commit()
    conn.close()

def update_protein_chembl_id(upac, chembl_id):
    conn = sqlite3.connect(PATH_DB)
    curs = conn.cursor()
    curs.execute('''
        UPDATE Protein
        SET ChemblId = ?
        WHERE UniProtAccNum = ?
        ''', (chembl_id, upac))
    conn.commit()
    conn.close()

def update_protein_gene_symbol(upac, gene_symbol):
    conn = sqlite3.connect(PATH_DB)
    curs = conn.cursor()
    curs.execute('''
        UPDATE Protein
        SET GeneSymbol = ?
        WHERE UniProtAccNum = ?
        ''', (gene_symbol, upac))
    conn.commit()
    conn.close()

def update_protein_process_function(upac, process_function):
    conn = sqlite3.connect(PATH_DB)
    curs = conn.cursor()
    curs.execute('''
        UPDATE Protein
        SET ProcessFunction = ?
        WHERE UniProtAccNum = ?
        ''', (process_function, upac))
    conn.commit()
    conn.close()

def update_protein_ions(upac, ions):
    conn = sqlite3.connect(PATH_DB)
    curs = conn.cursor()
    curs.execute('''
        UPDATE Protein
        SET Ions = ?
        WHERE UniProtAccNum = ?
        ''', (ions, upac))
    conn.commit()
    conn.close()

def update_protein_gating(upac, gating):
    conn = sqlite3.connect(PATH_DB)
    curs = conn.cursor()
    curs.execute('''
        UPDATE Protein SET Gating = ?
        WHERE UniProtAccNum = ?
        ''', (gating, upac))
    conn.commit()
    conn.close()

def update_protein_in_betse(upac, in_betse):
    conn = sqlite3.connect(PATH_DB)
    curs = conn.cursor()
    curs.execute('''
        UPDATE Protein
        SET InBETSE = ?
        WHERE UniProtAccNum = ?
        ''', (in_betse, upac))
    conn.commit()
    conn.close()

def update_protein_ion_channel_class_desc(upac, ion_channel_class_desc):
    conn = sqlite3.connect(PATH_DB)
    curs = conn.cursor()
    curs.execute('''
        UPDATE Protein
        SET IonChannelClassDesc = ?
        WHERE UniProtAccNum = ?
        ''', (ion_channel_class_desc, upac))
    conn.commit()
    conn.close()

def update_protein_ion_channel_sub_class(upac, ion_channel_sub_class):
    conn = sqlite3.connect(PATH_DB)
    curs = conn.cursor()
    curs.execute('''
        UPDATE Protein
        SET IonChannelSubClass = ?
        WHERE UniProtAccNum = ?
        ''', (ion_channel_sub_class, upac))
    conn.commit()
    conn.close()

def exists_protein(upac):
    conn = sqlite3.connect(PATH_DB)
    curs = conn.cursor()
    curs.execute('''
        SELECT EXISTS(
            SELECT 1
            FROM Protein
            WHERE UniProtAccNum=?
            LIMIT 1
        )
        ''', (upac,))
    row = curs.fetchone()
    if row[0] == 1:
        result = True
    else:
        result = False
    conn.close()
    return result

def lookup_protein(upac):
    conn = sqlite3.connect(PATH_DB)
    curs = conn.cursor()
    curs.execute('''
        SELECT UniProtAccNum, GeneSymbol, Name, ProcessFunction, Ions, Gating,
            InBETSE, IonChannelClassDesc, IonChannelSubClass, ChemblId
        FROM Protein
        WHERE UniProtAccNum=?
        ''', (upac,))
    resultset = curs.fetchone()
    if resultset == None:
        result = ['','','','','','','','','','']
    else:
        result = resultset
    conn.close()
    return result

def get_uniprots_by_gene_symbol(gene_symbol):
    conn = sqlite3.connect(PATH_DB)
    curs = conn.cursor()
    curs.execute('''
        SELECT UniProtAccNum
        FROM Protein
        WHERE GeneSymbol=?
        ''', (gene_symbol,))
    resultset = curs.fetchall()
    if resultset == None:
        result = []
    else:
        result = resultset
    conn.close()
    return resultset

def get_gene_symbol_by_uniprot(upac):
    conn = sqlite3.connect(PATH_DB)
    curs = conn.cursor()
    curs.execute('''
        SELECT GeneSymbol
        FROM Protein
        WHERE UniProtAccNum=?
        ''', (upac,))
    resultset = curs.fetchone()
    if resultset == None:
        result = ''
    else:
        result = resultset[0]
    conn.close()
    return result

def is_in_betse(upac):
    conn = sqlite3.connect(PATH_DB)
    curs = conn.cursor()
    curs.execute('''
        SELECT InBETSE
        FROM Protein
        WHERE UniProtAccNum=?
        ''', (upac,))
    resultset = curs.fetchone()
    if resultset == None:
        result = False
    else:
        if resultset[0] == 'Y':
            result = True
        else:
            result = False
    conn.close()
    return result

def get_all_protein_uniprots():
    conn = sqlite3.connect(PATH_DB)
    curs = conn.cursor()
    curs.execute('''
        SELECT UniProtAccNum
        FROM Protein
    ''')
    resultset = curs.fetchall()
    conn.close()
    return resultset

def dump_proteins():
    conn = sqlite3.connect(PATH_DB)
    curs = conn.cursor()
    with codecs.open('protein_dump.sql', 'w', 'utf-8') as f:
        for line in conn.iterdump():
            f.write('%s\n' % line)
    conn.close()

# Routines for Expression table

def add_expression(tissue_name, upac, expr_level, expr_level_qual='',
        expr_units='', assay_type='', dataset_name='', sourcedb_name=''):
    conn = sqlite3.connect(PATH_DB)
    curs = conn.cursor()
    curs.execute('''
        INSERT INTO Expression (TissueName, ProteinUniProtAccNum, ExprLevel,
            ExprLevelQual, ExprUnits, AssayType, DatasetName, SourceDBName)
        VALUES (?, ?, ?, ?, ?, ?, ?, ?)
        ''', (tissue_name, upac, expr_level, expr_level_qual, expr_units,
                assay_type, dataset_name, sourcedb_name))
    conn.commit()
    conn.close()

def exists_expression(upac):
    conn = sqlite3.connect(PATH_DB)
    curs = conn.cursor()
    curs.execute('''
        SELECT EXISTS(
            SELECT 1
            FROM Expression
            WHERE ProteinUniProtAccNum=?
            LIMIT 1
        )''', (upac,))
    row = curs.fetchone()
    if row[0] == 1:
        result = True
    else:
        result = False
    conn.close()
    return result

def lookup_expr(id):
    conn = sqlite3.connect(PATH_DB)
    curs = conn.cursor()
    curs.execute('''
        SELECT Id, TissueName, ProteinUniProtAccNum, ExprLevel, ExprLevelQual,
            ExprUnits, AssayType, DatasetName, SourceDBName
        FROM Expression
        WHERE Id=?
        ''', (id,))
    resultset = curs.fetchone()
    if resultset == None:
        result = ['', '', '', 0.0, '', '', '', '', '']
    else:
        result = resultset
    conn.close()
    return result
    
def get_expr_ids_by_uniprot(upac):
    conn = sqlite3.connect(PATH_DB)
    curs = conn.cursor()
    curs.execute('''
        SELECT Id
        FROM Expression
        WHERE ProteinUniProtAccNum=?
        ''', (upac,))
    resultset = curs.fetchall()
    conn.close()
    return resultset

def get_expr_for_tissue(tissue_name):
    prot_expr_list = []
    conn = sqlite3.connect(PATH_DB)
    curs = conn.cursor()
    for row in curs.execute('''
        SELECT ProteinUniProtAccNum, ExprLevel
        FROM Expression
        WHERE TissueName = ?
        ''',  (tissue_name,)):
            prot_expr_list.append((row[0], row[1]))
    conn.commit()
    conn.close()
    return prot_expr_list

def get_expr_levels(upac, tissue):
    conn = sqlite3.connect(PATH_DB)
    curs = conn.cursor()
    curs.execute('''
        SELECT ExprLevel
        FROM Expression
        WHERE ProteinUniProtAccNum=? AND TissueName=?
        ''', (upac, tissue))
    resultset = curs.fetchall()
    conn.close()
    return resultset

def get_expr_level_quals(upac, tissue):
    conn = sqlite3.connect(PATH_DB)
    curs = conn.cursor()
    curs.execute('''
        SELECT ExprLevelQual
        FROM Expression
        WHERE ProteinUniProtAccNum=? AND TissueName=?
        AND ExprLevelQual!=""
        ''', (upac, tissue))
    resultset = curs.fetchall()
    conn.close()
    return resultset

def get_expr_level_for_dataset(upac, tissue, dataset):
    conn = sqlite3.connect(PATH_DB)
    curs = conn.cursor()
    curs.execute('''
        SELECT ExprLevel
        FROM Expression
        WHERE ProteinUniProtAccNum=? AND TissueName=? AND DatasetName=?
        ''',
        (upac, tissue, dataset))
    resultset = curs.fetchall()
    conn.close()
    return resultset

def get_all_expression():
    conn = sqlite3.connect(PATH_DB)
    curs = conn.cursor()
    curs.execute('''
        SELECT Id
        FROM Expression
        ''')
    resultset = curs.fetchall()
    conn.close()
    return resultset

def exists_expr_threshold(tissue, upac, threshold, dataset):
    conn = sqlite3.connect(PATH_DB)
    curs = conn.cursor()
    curs.execute('''
        SELECT EXISTS(
            SELECT 1
            FROM Expression
            WHERE (
                TissueName=? AND ProteinUniProtAccNum=? AND DatasetName=?
                AND ExprLevel>=?
            )
            LIMIT 1
        )
        ''', (tissue, upac, dataset, threshold))
    row = curs.fetchone()
    if row[0] == 1:
        result = True
    else:
        result = False
    conn.close()
    return result

# Routines for DBGene table

def add_dbgene(externaldb_name, gbac, probe_id):
    conn = sqlite3.connect(PATH_DB)
    curs = conn.cursor()
    curs.execute('''
        INSERT INTO DBGene
        (ExternalDBName, GenbankAccNum, ProbeID)
        VALUES (?, ?, ?)
        ''', (externaldb_name, gbac, probe_id))
    conn.commit()
    conn.close()

# Routines for GenbankUniprot table

def add_genbank_uniprot(gbac, upac):
    conn = sqlite3.connect(PATH_DB)
    curs = conn.cursor()
    curs.execute('''
        INSERT INTO GenbankUniprot
        (GenbankAccNum, UniProtAccNum)
        VALUES (?, ?)
        ''', (gbac, upac))
    conn.commit()
    conn.close()    

def exists_genbank_uniprot(gbac, upac):
    conn = sqlite3.connect(PATH_DB)
    curs = conn.cursor()
    curs.execute('''
        SELECT EXISTS(
            SELECT 1
            FROM GenbankUniprot
            WHERE GenbankAccNum=? AND UniProtAccNum=?
            LIMIT 1
        )
        ''', (gbac, upac))
    row = curs.fetchone()
    if row[0] == 1:
        result = True
    else:
        result = False
    conn.close()
    return result

# Routines for Compound table

def add_compound(name, smiles='', inchi='', chembl_id='', synonyms='',
        approval_status='', first_approval_year='', sourcedb_name=''):
    conn = sqlite3.connect(PATH_DB)
    curs = conn.cursor()
    curs.execute('''
        INSERT INTO Compound
        (SMILES, InChI, Name, ChemblId, Synonyms, ApprovalStatus,
            FirstApprovalYear, SourceDBName)
        VALUES (?, ?, ?, ?, ?, ?, ?, ?)
        ''', (smiles, inchi, name, chembl_id, synonyms, approval_status,
                first_approval_year, sourcedb_name))
    conn.commit()
    conn.close()    

def exists_compound_by_chembl_id(chembl_id):
    conn = sqlite3.connect(PATH_DB)
    curs = conn.cursor()
    expr_list = []
    curs.execute('''
        SELECT EXISTS(
            SELECT 1
            FROM Compound
            WHERE ChemblId=?
            LIMIT 1
        )
        ''', (chembl_id,))
    row = curs.fetchone()
    if row[0] == 1:
        result = True
    else:
        result = False
    conn.close()
    return result

def lookup_compound(id):
    conn = sqlite3.connect(PATH_DB)
    curs = conn.cursor()
    curs.execute('''
        SELECT Id, SMILES, InChI, Name, ChemblId, Synonyms, ApprovalStatus,
            FirstApprovalYear, SourceDBName
        FROM Compound WHERE Id=?
        ''', (id,))
    resultset = curs.fetchone()
    if resultset == None:
        result = ['', '', '', '', '', '', '', '', '']
    else:
        result = resultset
    conn.close()
    return result

def get_compound_id_by_chembl_id(chembl_id):
    conn = sqlite3.connect(PATH_DB)
    curs = conn.cursor()
    curs.execute('''
        SELECT Id
        FROM Compound
        WHERE ChemblId=?
        ''', (chembl_id,))
    resultset = curs.fetchone()
    if resultset == None:
        result = ''
    else:
        result = resultset[0]
    conn.close()
    return result

def get_compound_id_by_compound_name(compound_name):
    conn = sqlite3.connect(PATH_DB)
    curs = conn.cursor()
    curs.execute('''
        SELECT Id
        FROM Compound
        WHERE Name=?
        ''', (compound_name,))
    resultset = curs.fetchone()
    if resultset == None:
        result = ''
    else:
        result = resultset[0]
    conn.close()
    return result

# Routines for Interaction table

def add_interaction(target_upac, compound_id, action_type='',
        action_desc='', strength=0, strength_units='', assay_type='',
        chembl_id='', sourcedb_name=''):
    conn = sqlite3.connect(PATH_DB)
    curs = conn.cursor()
    curs.execute('''
        INSERT INTO Interaction
        (TargetUniProtAccNum, CompoundId, ActionType, ActionDesc, Strength,
            StrengthUnits, AssayType, ChemblId, SourceDBName)
        VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?)
        ''', (target_upac, compound_id, action_type, action_desc, strength,
                strength_units, assay_type, chembl_id, sourcedb_name))
    conn.commit()
    conn.close()

def exists_interaction(upac, compound_id):
    conn = sqlite3.connect(PATH_DB)
    curs = conn.cursor()
    curs.execute('''
        SELECT EXISTS(
            SELECT 1
            FROM Interaction
            WHERE TargetUniProtAccNum=? AND CompoundId=?
            LIMIT 1
        )
        ''', (upac, compound_id))
    row = curs.fetchone()
    if row[0] == 1:
        result = True
    else:
        result = False
    conn.close()
    return result

def exists_any_interaction_for_protein(upac):
    conn = sqlite3.connect(PATH_DB)
    curs = conn.cursor()
    curs.execute('''
        SELECT EXISTS(
            SELECT 1
            FROM Interaction
            WHERE TargetUniProtAccNum=?
            LIMIT 1
        )
        ''', (upac,))
    row = curs.fetchone()
    if row[0] == 1:
        result = True
    else:
        result = False
    conn.close()
    return result

def get_interaction_ids_by_uniprot(upac):
    conn = sqlite3.connect(PATH_DB)
    curs = conn.cursor()
    curs.execute('''
        SELECT Id
        FROM Interaction
        WHERE TargetUniProtAccNum=?
        ''', (upac,))
    resultset = curs.fetchall()
    conn.close()
    return resultset

def lookup_interaction(id):
    conn = sqlite3.connect(PATH_DB)
    curs = conn.cursor()
    curs.execute('''
        SELECT Id, TargetUniProtAccNum, CompoundId, ActionType, ActionDesc,
            Strength, StrengthUnits, AssayType, ChemblId, SourceDBName
        FROM Interaction
        WHERE Id=?
        ''', (id,))
    resultset = curs.fetchone()
    if resultset == None:
        result = ['', '', '', '', '', '', '', '', '', '']
    else:
        result = resultset[0]
    conn.close()
    return resultset

# Routines for DBTissue table

def add_dbtissue(externaldb_name, tissue_name, db_equivalent_tissue_name):
    conn = sqlite3.connect(PATH_DB)
    curs = conn.cursor()
    curs.execute('''
        INSERT INTO DBTissue
        (ExternalDBName, TissueName, DBEquivalentTissueName)
        VALUES (?, ?, ?)
        ''', (externaldb_name, tissue_name, db_equivalent_tissue_name))
    conn.commit()
    conn.close()
    
def get_bto_tissue_in_dbtissue(externaldb_name, db_equivalent_tissue_name):
    conn = sqlite3.connect(PATH_DB)
    curs = conn.cursor()
    curs.execute('''
        SELECT TissueName
        FROM DBTissue
        WHERE ExternalDBName=? AND DBEquivalentTissueName=?
        ''', (externaldb_name, db_equivalent_tissue_name))
    resultset = curs.fetchone()
    if resultset == None:
        result = ''
    else:
        result = resultset[0]
    conn.close()
    return result

# Routines for Specificity table

def add_specificity(tissue_name, upac, specificity_score):
    conn = sqlite3.connect(PATH_DB)
    curs = conn.cursor()
    curs.execute('''
        INSERT INTO Specificity
        (TissueName, UniProtAccNum, SpecificityScore)
        VALUES (?, ?, ?)
        ''', (tissue_name, upac, specificity_score))
    conn.commit()
    conn.close()

# Routines for multiple tables

def dump_database():
    conn = sqlite3.connect(PATH_DB)
    with codecs.open('levindb_dump.sql', 'w', 'utf-8') as f:
        for line in conn.iterdump():
            f.write('%s\n' % line)
    conn.close()

def print_db_stats():
    conn = sqlite3.connect(PATH_DB)
    curs = conn.cursor()
    curs.execute("SELECT * FROM sqlite_master WHERE type='table'")
    table_name_list = []
    for table_info in curs:
        table_name = table_info[1]
        table_name_list.append(table_name)
    conn.close()
    for table_name in table_name_list:
        conn = sqlite3.connect(PATH_DB)
        curs = conn.cursor()
        curs.execute('SELECT COUNT(*) FROM %s' % (table_name,))
        row_count = curs.fetchone()[0]
        print table_name, row_count
        conn.close()
    useful_data_count = 0
    upac_record_list = get_all_protein_uniprots()
    for upac_record in upac_record_list:
        upac = upac_record[0]
        interaction_id_list = get_interaction_ids_by_uniprot(upac)
        expression_id_list = get_expr_ids_by_uniprot(upac)
        if len(interaction_id_list) > 0 and len(expression_id_list) > 0:
            useful_data_count += 1
    print 'Protein with interaction and expression %d' % useful_data_count
    upac_record_list = get_all_protein_uniprots()
    in_betse_count = 0
    in_betse_with_expr_count = 0
    list_in_betse_with_expr = []
    list_in_betse_no_expr = []
    for upac_record in upac_record_list:
        upac = upac_record[0]
        protein_record = lookup_protein(upac)
        in_betse = protein_record[6]
        if in_betse == 'Y':
            in_betse_count += 1
            if exists_expression(upac):
                in_betse_with_expr_count += 1
                list_in_betse_with_expr.append(get_gene_symbol_by_uniprot(upac))
            else:
                list_in_betse_no_expr.append(get_gene_symbol_by_uniprot(upac))
    print 'Number of proteins in betse = %d' % in_betse_count
    print '...also with with expression = %d' % in_betse_with_expr_count

# DB management functions

def cleanup_externaldb():
    conn = sqlite3.connect(PATH_DB)
    curs = conn.cursor()
    curs.execute('''
        DELETE FROM ExternalDB
        ''')
    conn.commit()
    conn.close()

def setup_externaldb():
    add_externaldb(EXTDB_BIOGPS, 'http://biogps.org')
    add_externaldb(EXTDB_CHEMBL, 'https://www.ebi.ac.uk/chembl/')
    add_externaldb(EXTDB_BRENDA, 'http://www.brenda-enzymes.org')
    add_externaldb(EXTDB_DRUGBANK, 'https://www.drugbank.ca')
    add_externaldb(EXTDB_ZINC15, 'http://zinc15.docking.org')
    add_externaldb(EXTDB_HMDB, 'http://www.hmdb.ca')
    add_externaldb(EXTDB_CHANNELPEDIA, 'http://channelpedia.epfl.ch')

def cleanup_ontology():
    conn = sqlite3.connect(PATH_DB)
    curs = conn.cursor()
    curs.execute('''
        DELETE FROM Tissue
        ''')
    conn.commit()
    conn.close()

def setup_ontology():
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
            add_tissue(tissue_name, bto_id, '')

def cleanup_genbank_uniprot():
    conn = sqlite3.connect(PATH_DB)
    curs = conn.cursor()
    curs.execute('''
        DELETE FROM GenbankUniprot
        ''')
    conn.commit()
    conn.close()

def setup_genbank_uniprot():
    gu_file = open(PATH_UNIPROT_GENBANK, 'rU')
    line_num = 0
    for line in gu_file:
        if line_num != 0:
            line_split = line.split('\t')
            if len(line_split) == 2:
                gbac = line_split[0].strip()
                upac = line_split[1].strip()
                if upac != '':
                    if not exists_genbank_uniprot(gbac, upac):
                        add_genbank_uniprot(gbac, upac)
        line_num += 1
    gu_file.close()

def cleanup_proteins():
    conn = sqlite3.connect(PATH_DB)
    curs = conn.cursor()
    curs.execute('''
        DELETE FROM Protein
        ''')
    conn.commit()
    conn.close()

def setup_ion_channels():
    prot_file = open(PATH_GO_ION_CHANNELS, 'rU')
    prot_reader = csv.reader(prot_file)
    prot_reader.next()
    for row in prot_reader:
        protein_name = row[2]
        process_function = row[3]
        upac = row[0]
        if upac == '':
            continue
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
        if not exists_protein(upac):
            add_protein(upac, gene_symbol, protein_name, process_function,
                    ions, gating, in_betse, ion_channel_subclass)
        else:
            old_protein_record = lookup_protein(upac)
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
            update_protein_gene_symbol(upac, gene_symbol)
            update_protein_name(upac, protein_name)
            update_protein_ions(upac, ions)
            update_protein_gating(upac, gating)
            if old_in_betse == '':
                update_protein_in_betse(upac, in_betse)
            update_protein_ion_channel_sub_class(upac, ion_channel_subclass)
    prot_file.close()

def cleanup_channel_classes():
    conn = sqlite3.connect(PATH_DB)
    curs = conn.cursor()
    curs.execute('''
        DELETE FROM ChannelSuperClass
        ''')
    conn.commit()
    conn.close()
    conn = sqlite3.connect(PATH_DB)
    curs = conn.cursor()
    curs.execute('''
        DELETE FROM ChannelClass
        ''')
    conn.commit()
    conn.close()
    conn = sqlite3.connect(PATH_DB)
    curs = conn.cursor()
    curs.execute('''
        DELETE FROM ChannelSubClass
        ''')
    conn.commit()
    conn.close()

def setup_channel_classes():
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
        if not exists_channel_super_class(channel_superclass):
            add_channel_super_class(channel_superclass)
        if not exists_channel_class(channel_class):
            add_channel_class(channel_class, channel_superclass)
        if not exists_channel_sub_class(channel_subclass):
            if channel_subclass in channelpedia_dict:
                channelpedia_info, channelpedia_url = \
                        channelpedia_dict[channel_subclass]
            else:
                channelpedia_info = ''
                channelpedia_url = ''
            add_channel_sub_class(channel_subclass, channel_class,
                    channelpedia_text=channelpedia_info.decode('utf_8').encode('utf_8'),
                    channelpedia_url=channelpedia_url.decode('utf_8').encode('utf_8'))
        upac_record_list = get_uniprots_by_gene_symbol(gene_symbol)
        for upac_record in upac_record_list:
            upac = upac_record[0]
            update_protein_sub_class(upac, channel_subclass)
    class_file.close()

def cleanup_dbtissue():
    conn = sqlite3.connect(PATH_DB)
    curs = conn.cursor()
    curs.execute('''
        DELETE FROM DBTissue
        ''')
    conn.commit()
    conn.close()

def setup_dbtissue():
    setup_biogps()
    setup_hpa()

def setup_biogps():
    translation_file = open(PATH_BIOGPS_TRANSLATION_TABLE, 'rU')
    translation_reader = csv.reader(translation_file)
    translation_reader.next()
    for row in translation_reader:
        biogps_tissue = row[0]
        bto_tissue = row[1]
        add_dbtissue(EXTDB_BIOGPS, bto_tissue, biogps_tissue)
    translation_file.close()

def setup_hpa():
    translation_file = open(PATH_HPA_TRANSLATION_TABLE, 'rU')
    translation_reader = csv.reader(translation_file)
    translation_reader.next()
    for row in translation_reader:
        hpa_tissue = row[0]
        bto_tissue = row[1]
        add_dbtissue(EXTDB_HPA, bto_tissue, hpa_tissue)
    translation_file.close()

def cleanup_expression():
    conn = sqlite3.connect(PATH_DB)
    curs = conn.cursor()
    curs.execute('''
        DELETE FROM Expression
        ''')
    conn.commit()
    conn.close()

def read_biogps_annot():
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

def process_biogps_data():
    annot_dict = read_biogps_annot()
    microarray_file = open(PATH_BIOGPS_GCRMA, 'rU')
    microarray_reader = csv.reader(microarray_file)
    header = microarray_reader.next()
    tissue_list = header[1:]
    for row in microarray_reader:
        probeset_id = row[0]
        if not probeset_id in annot_dict:
            continue
        gene_symbol = annot_dict[probeset_id]
        upac_resultset = get_uniprots_by_gene_symbol(gene_symbol.strip())
        if len(upac_resultset) == 0:
            continue
        for i in xrange(0, len(tissue_list)):
            expr_level = float(row[i + 1])
            tissue_name = get_bto_tissue_in_dbtissue(EXTDB_BIOGPS, 
                    tissue_list[i])
            for upac_result in upac_resultset:
                upac = upac_result[0]
                add_expression(tissue_name, upac, expr_level, 
                        assay_type=ASSAY_MICROARRAY,
                        dataset_name=DATASET_BIOGPS, 
                        sourcedb_name=EXTDB_BIOGPS)
    microarray_file.close()

def read_hpa_map():
    ensembl_to_uniprot = {}
    map_file = open(PATH_UNIPROT_ENSEMBL_IDMAPPING)
    for line in map_file:
        row = line.strip().split('\t')
        upac = row[0]
        ensg = row[2]
        ensembl_to_uniprot[ensg] = upac
    map_file.close()
    return ensembl_to_uniprot

def process_hpa_rna_data():
    hpa_file = open(PATH_HPA_RNA_TISSUE)
    hpa_file.next()
    for line in hpa_file:
        row = line.strip().split('\t')
        if len(row) >= 5 and row[0].strip() != '':
            gene_symbol = row[1]
            upac_record_list = get_uniprots_by_gene_symbol(gene_symbol)
            if len(upac_record_list) != 0:
                for upac_record in upac_record_list:
                    upac = upac_record[0]
                    gene_name = row[1]
                    tissue_name = row[2]
                    tissue_bto = get_bto_tissue_in_dbtissue(EXTDB_HPA, 
                            tissue_name)
                    expr_value = float(row[3])
                    add_expression(tissue_bto, upac, expr_value,
                            assay_type=ASSAY_RNASEQ, dataset_name=DATASET_HPA,
                            sourcedb_name=EXTDB_HPA)
    hpa_file.close()

def process_hpa_cancer_data():
    hpa_file = open(PATH_HPA_PATHOLOGY)
    hpa_file.next()
    expr_qual_choices = ['High', 'Medium', 'Low', 'Not detected']
    for line in hpa_file:
        row = line.strip().split('\t')
        if len(row) >= 6:
            gene_symbol = row[1]
            upac_record_list = get_uniprots_by_gene_symbol(gene_symbol)
            if len(upac_record_list) != 0:
                for upac_record in upac_record_list:
                    upac = upac_record[0]
                    tissue_name = row[2]
                    tissue_bto = get_bto_tissue_in_dbtissue(EXTDB_HPA,
                            tissue_name)
                    expr_value = 0.0
                    data_is_okay = False
                    try:
                         patient_counts = map(int, row[3:7])
                         data_is_okay = True
                    except:
                        pass
                    if data_is_okay:
                        patient_count_max = max(patient_counts)
                        patient_count_max_idx = [idx for idx, value in enumerate(patient_counts) if value==patient_count_max][-1]
                        expr_value_qual = \
                                expr_qual_choices[patient_count_max_idx]
                        add_expression(tissue_bto, upac, expr_value,
                            expr_level_qual=expr_value_qual,
                            assay_type=ASSAY_IMMUNO, dataset_name=DATASET_HPA,
                            sourcedb_name=EXTDB_HPA)
    hpa_file.close()

def process_hpa_normal_data():
    hpa_file = open(PATH_HPA_PATHOLOGY)
    hpa_file.next()
    for line in hpa_file:
        row = line.strip().split('\t')
        if len(row) >= 6:
            gene_symbol = row[1]
            upac_record_list = get_uniprots_by_gene_symbol(gene_symbol)
            if len(upac_record_list) == 6:
                for upac_record in upac_record_list:
                    upac = upac_record[0]
                    tissue_name = row[2]
                    tissue_bto = get_bto_tissue_in_dbtissue(EXTDB_HPA, 
                            tissue_name)
                    expr_value = 0.0
                    expr_value_qual = row[4]
                    if expr_value_qual in ['Not detected', 'Low', 'Medium', 
                            'High']:
                        add_expression(tissue_bto, upac, expr_value,
                            expr_level_qual=expr_value_qual,
                            assay_type=ASSAY_IMMUNO, dataset_name=DATASET_HPA,
                            sourcedb_name=EXTDB_HPA)
    hpa_file.close()

def cleanup_compounds():
    conn = sqlite3.connect(PATH_DB)
    curs = conn.cursor()
    curs.execute('''DELETE FROM Compound''')
    conn.commit()
    conn.close()

def cleanup_interactions():
    conn = sqlite3.connect(PATH_DB)
    curs = conn.cursor()
    curs.execute('''DELETE FROM Interaction''')
    conn.commit()
    conn.close()

def input_chembl_assays():
    input_chembl_compounds()
    input_chembl_drugs()

def input_chembl_compounds():
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
        if not exists_compound_by_chembl_id(compound_chembl_id):
            add_compound(compound_name, '', '', chembl_id=compound_chembl_id,
                    synonyms=iupac_name, sourcedb_name=EXTDB_CHEMBL)
        if not exists_protein(upac):
            add_protein(upac, name=target_name, process_function=target_type,
            chembl_id=target_chembl_id)
        compound_id = get_compound_id_by_chembl_id(compound_chembl_id)
        add_interaction(upac, compound_id, action_type='',
                strength=assay_value.decode('utf_8').encode('utf_8'),
                strength_units=assay_units.decode('utf_8').encode('utf_8'),
                assay_type=assay_standard_type.decode('utf_8').encode('utf_8'),
                chembl_id=assay_chembl_id, sourcedb_name=EXTDB_CHEMBL)
    compound_file.close()

def input_chembl_drugs():
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
        if not exists_compound_by_chembl_id(active_compound_chembl_id):
            add_compound('', '', active_compound_name,
                    chembl_id=active_compound_chembl_id,
                    synonyms=dosed_compound_name, sourcedb_name=EXTDB_CHEMBL)
        if not exists_protein(upac):
            add_protein(upac, name=target_name, process_function=target_type,
                    chembl_id=target_chembl_id)
        compound_id = get_compound_id_by_chembl_id(active_compound_chembl_id)
        add_interaction(upac, compound_id,
                action_type=drug_mechanism.decode('utf_8').encode('utf_8'),
                strength=assay_value.decode('utf_8').encode('utf_8'),
                strength_units=assay_units.decode('utf_8').encode('utf_8'),
                assay_type=assay_standard_type.decode('utf_8').encode('utf_8'),
                chembl_id=assay_chembl_id, sourcedb_name=EXTDB_CHEMBL)
    drug_file.close()

def cleanup_specificity():
    conn = sqlite3.connect(PATH_DB)
    curs = conn.cursor()
    curs.execute('''
        DELETE FROM Specificity
        ''')
    conn.commit()
    conn.close()

def calculate_specificity():
    UBIQUITY_PARAM = 0.3
    expr_table = {}
    ubiquity_table = {}
    specificity_table = {}
    all_upac_list = []
    all_tissue_name_list = []
    upac_record_list = get_all_protein_uniprots()
    for upac_record in upac_record_list:
        upac = upac_record[0]
        all_upac_list.append(upac)
    tissue_name_record_list = get_all_tissue_names()
    for tissue_name_record in tissue_name_record_list:
        tissue_name = tissue_name_record[0]
        all_tissue_name_list.append(tissue_name)
    for tissue_name in all_tissue_name_list:
        for upac in all_upac_list:
            expr_level_qual_record_list = get_expr_level_quals(upac, 
                    tissue_name)
            expr_level_qual_list = []
            for expr_level_qual_record in expr_level_qual_record_list:
                expr_level_qual = expr_level_qual_record[0]
                expr_level_qual_list.append(expr_level_qual)
            if len(expr_level_qual_list) != 0:
                if 'High' in expr_level_qual_list:
                    expr_level_qual = 'High'
                elif 'Medium' in expr_level_qual_list:
                    expr_level_qual = 'Medium'
                elif 'Low' in expr_level_qual_list:
                    expr_level_qual = 'Low'
                else:
                    expr_level_qual = 'Not detected'
                expr_table[(tissue_name, upac)] = expr_level_qual
        #print 'Checked %s' % (tissue_name,)
    for upac in all_upac_list:
        ubiquity_score = 0
        for tissue_name in all_tissue_name_list:
            if (tissue_name, upac) in expr_table:
                expr_level_qual = expr_table[(tissue_name, upac)]
                if expr_level_qual in ['High', 'Medium']:
                    ubiquity_score += 1
        ubiquity_table[upac] = ubiquity_score
    ubiquity_threshold = round(UBIQUITY_PARAM * len(all_tissue_name_list))
    for upac in all_upac_list:
        ubiquity_score = ubiquity_table[upac]
        for tissue_name in all_tissue_name_list:
            if (tissue_name, upac) in expr_table:
                expr_level_qual = expr_table[(tissue_name, upac)]
                if expr_level_qual == 'High' and \
                        ubiquity_score <= ubiquity_threshold:
                    specificity_score = '1'
                else:
                    specificity_score = '0'
                specificity_table[(tissue_name, upac)] = specificity_score
                add_specificity(tissue_name, upac, specificity_score)

def calculate_specificity_jp():
    all_upac_list = []
    upac_record_list = get_all_protein_uniprots()
    for upac_record in upac_record_list:
        upac = upac_record[0]
        all_upac_list.append(upac)
    tissue_name_record_list = get_all_tissue_names()
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
            expr_level_record_list = get_expr_level_for_dataset(upac, 
                    tissue_name, EXTDB_HPA)
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
            add_specificity(tissue_name, upac, specificity_score_tissue)
                
