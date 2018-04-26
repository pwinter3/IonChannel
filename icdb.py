import codecs
import sqlite3


class IonChannelDatabase(object):

    def __init__(self, path_db):
        self.path_db = path_db
        self.conn = sqlite3.connect(path_db)

    def get_conn(self):
        return self.conn


    # Routines for ChannelSuperClass table

    def create_channel_super_class_table(self):
        conn = self.get_conn()
        curs = conn.cursor()
        curs.execute("""
            DROP TABLE IF EXISTS ChannelSuperClass
            """)
        curs.execute("""
            CREATE TABLE ChannelSuperClass(
                Name TEXT PRIMARY KEY NOT NULL)
            """)
        conn.commit()

    def add_channel_super_class(self, name):
        conn = self.get_conn()
        curs = conn.cursor()
        curs.execute("""
            INSERT INTO ChannelSuperClass(Name)
            VALUES (?)
            """, (name,))
        conn.commit()

    def exists_channel_super_class(self, name):
        conn = self.get_conn()
        curs = conn.cursor()
        curs.execute("""
            SELECT EXISTS(
                SELECT 1
                FROM ChannelSuperClass
                WHERE Name = ?
                LIMIT 1)
            """, (name,))
        row = curs.fetchone()
        if row[0] == 1:
            result = True
        else:
            result = False
        return result

    def delete_channel_super_classes(self):
        conn = self.get_conn()
        curs = conn.cursor()
        curs.execute("""
            DELETE FROM ChannelSuperClass
            """)
        conn.commit()


    # Routines for ChannelClass table

    def create_channel_class_table(self):
        conn = self.get_conn()
        curs = conn.cursor()
        curs.execute("""
            DROP TABLE IF EXISTS ChannelClass
            """)
        curs.execute("""
            CREATE TABLE ChannelClass(
                Name TEXT PRIMARY KEY NOT NULL,
                SuperClass TEXT NOT NULL)
            """)
        conn.commit()

    def add_channel_class(self, name, super_class=''):
        conn = self.get_conn()
        curs = conn.cursor()
        curs.execute("""
            INSERT INTO ChannelClass(Name, SuperClass)
            VALUES (?, ?)
            """, (name, super_class))
        conn.commit()

    def exists_channel_class(self, name):
        conn = self.get_conn()
        curs = conn.cursor()
        curs.execute("""
            SELECT EXISTS(
                SELECT 1
                FROM ChannelClass
                WHERE Name = ?
                LIMIT 1)
            """, (name,))
        row = curs.fetchone()
        if row[0] == 1:
            result = True
        else:
            result = False
        return result

    def delete_channel_classes(self):
        conn = self.get_conn()
        curs = conn.cursor()
        curs.execute("""
            DELETE FROM ChannelClass
            """)
        conn.commit()


    # Routines for ChannelSubClass table

    def create_channel_sub_class_table(self):
        conn = self.get_conn()
        curs = conn.cursor()
        curs.execute("""
            DROP TABLE IF EXISTS ChannelSubClass
            """)
        curs.execute("""
            CREATE TABLE ChannelSubClass(
                Name TEXT PRIMARY KEY NOT NULL,
                Class TEXT NOT NULL,
                ChannelpediaText TEXT NOT NULL,
                ChannelpediaURL TEXT NOT NULL)
            """)
        conn.commit()

    def add_channel_sub_class(
            self, name, channel_class='', channelpedia_text='',
            channelpedia_url=''):
        conn = self.get_conn()
        conn.text_factory = str
        curs = conn.cursor()
        curs.execute("""
            INSERT INTO ChannelSubClass(
                Name, Class, ChannelpediaText, ChannelpediaURL)
            VALUES (?, ?, ?, ?)
            """, (name, channel_class, channelpedia_text, channelpedia_url))
        conn.commit()

    def exists_channel_sub_class(self, name):
        conn = self.get_conn()
        curs = conn.cursor()
        expr_list = []
        curs.execute("""
            SELECT EXISTS(
                SELECT 1
                FROM ChannelSubClass
                WHERE Name = ?
                LIMIT 1)
            """, (name,))
        row = curs.fetchone()
        if row[0] == 1:
            result = True
        else:
            result = False
        return result

    def delete_channel_sub_classes(self):
        conn = self.get_conn()
        curs = conn.cursor()
        curs.execute("""
            DELETE FROM ChannelSubClass
            """)
        conn.commit()


    # Routines for ExternalDB table

    def create_external_db_table(self):
        conn = self.get_conn()
        curs = conn.cursor()
        curs.execute("""
            DROP TABLE IF EXISTS ExternalDB
            """)
        curs.execute("""
            CREATE TABLE ExternalDB(
                Name TEXT PRIMARY KEY NOT NULL,
                URL TEXT NOT NULL)
            """)
        conn.commit()

    def add_external_db(self, name, url=''):
        conn = self.get_conn()
        curs = conn.cursor()
        curs.execute("""
            INSERT INTO ExternalDB(Name, URL)
            VALUES (?, ?)
            """, (name, url))
        conn.commit()

    def delete_external_dbs(self):
        conn = sqlite3.connect(PATH_DB)
        curs = conn.cursor()
        curs.execute("""
            DELETE FROM ExternalDB
            """)
        conn.commit()


    # Routines for Tissue table

    def create_tissue_table(self):
        conn = self.get_conn()
        curs = conn.cursor()
        curs.execute("""
            DROP TABLE IF EXISTS Tissue
            """)
        curs.execute("""
            CREATE TABLE Tissue(
                Name TEXT PRIMARY KEY NOT NULL,
                BTOId TEXT UNIQUE NOT NULL,
                ParentTissueName TEXT NOT NULL)
            """)
        conn.commit()

    def add_tissue(self, name, bto_id='', parent_tissue_name=''):
        conn = self.get_conn()
        curs = conn.cursor()
        curs.execute("""
            INSERT INTO Tissue(Name, BTOId, ParentTissueName)
            VALUES (?, ?, ?)
            """, (name, bto_id, parent_tissue_name))
        conn.commit()

    def get_tissue_names(self):
        conn = self.get_conn()
        curs = conn.cursor()
        curs.execute("""
            SELECT Name
            FROM Tissue
            """)
        resultset = curs.fetchall()
        conn = self.get_conn()
        return resultset

    def dump_tissue_table(self):
        conn = self.get_conn()
        curs = conn.cursor()
        with codecs.open('tissue_dump.sql', 'w', 'utf-8') as f:
            for line in conn.iterdump():
                f.write('%s\n' % line)

    def delete_tissues(self):
        conn = self.get_conn()
        curs = conn.cursor()
        curs.execute("""
            DELETE FROM Tissue
            """)
        conn.commit()


    # Routines for Protein table

    def create_protein_table(self):
        conn = self.get_conn()
        curs = conn.cursor()
        curs.execute("""
            DROP TABLE IF EXISTS Protein
            """)
        curs.execute("""
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
                ChemblId TEXT NOT NULL)
            """)
        conn.commit()

    def add_protein(
            self, upac, gene_symbol='', name='', process_function='', ions='',
            gating='', in_betse='', ion_channel_class_desc='',
            ion_channel_sub_class='', chembl_id=''):
        conn = self.get_conn()
        curs = conn.cursor()
        curs.execute("""
            INSERT INTO Protein(
                UniProtAccNum, GeneSymbol, Name, ProcessFunction, Ions, Gating,
                InBETSE, IonChannelClassDesc, IonChannelSubClass, ChemblId)
            VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
            """, (upac, gene_symbol, name, process_function, ions, gating,
            in_betse, ion_channel_class_desc, ion_channel_sub_class, chembl_id))
        conn.commit()

    def update_protein_sub_class(self, upac, sub_class):
        conn = self.get_conn()
        curs = conn.cursor()
        curs.execute("""
            UPDATE Protein
            SET IonChannelSubClass = ?
            WHERE UniProtAccNum = ?
            """, (sub_class, upac))
        conn.commit()

    def update_protein_name(self, upac, name):
        conn = self.get_conn()
        curs = conn.cursor()
        curs.execute("""
            UPDATE Protein
            SET Name = ?
            WHERE UniProtAccNum = ?
            """, (name, upac))
        conn.commit()

    def update_protein_chembl_id(self, upac, chembl_id):
        conn = self.get_conn()
        curs = conn.cursor()
        curs.execute("""
            UPDATE Protein
            SET ChemblId = ?
            WHERE UniProtAccNum = ?
            """, (chembl_id, upac))
        conn.commit()

    def update_protein_gene_symbol(self, upac, gene_symbol):
        conn = self.get_conn()
        curs = conn.cursor()
        curs.execute("""
            UPDATE Protein
            SET GeneSymbol = ?
            WHERE UniProtAccNum = ?
            """, (gene_symbol, upac))
        conn.commit()

    def update_protein_process_function(self, upac, process_function):
        conn = self.get_conn()
        curs = conn.cursor()
        curs.execute("""
            UPDATE Protein
            SET ProcessFunction = ?
            WHERE UniProtAccNum = ?
            """, (process_function, upac))
        conn.commit()

    def update_protein_ions(self, upac, ions):
        conn = self.get_conn()
        curs = conn.cursor()
        curs.execute("""
            UPDATE Protein
            SET Ions = ?
            WHERE UniProtAccNum = ?
            """, (ions, upac))
        conn.commit()

    def update_protein_gating(self, upac, gating):
        conn = self.get_conn()
        curs = conn.cursor()
        curs.execute("""
            UPDATE Protein SET Gating = ?
            WHERE UniProtAccNum = ?
            """, (gating, upac))
        conn.commit()

    def update_protein_in_betse(self, upac, in_betse):
        conn = self.get_conn()
        curs = conn.cursor()
        curs.execute("""
            UPDATE Protein
            SET InBETSE = ?
            WHERE UniProtAccNum = ?
            """, (in_betse, upac))
        conn.commit()

    def update_protein_ion_channel_class_desc(self, upac,
            ion_channel_class_desc):
        conn = self.get_conn()
        curs = conn.cursor()
        curs.execute("""
            UPDATE Protein
            SET IonChannelClassDesc = ?
            WHERE UniProtAccNum = ?
            """, (ion_channel_class_desc, upac))
        conn.commit()

    def update_protein_ion_channel_sub_class(
            self, upac, ion_channel_sub_class):
        conn = self.get_conn()
        curs = conn.cursor()
        curs.execute("""
            UPDATE Protein
            SET IonChannelSubClass = ?
            WHERE UniProtAccNum = ?
            """, (ion_channel_sub_class, upac))
        conn.commit()

    def exists_protein(self, upac):
        conn = self.get_conn()
        curs = conn.cursor()
        curs.execute("""
            SELECT EXISTS(
                SELECT 1
                FROM Protein
                WHERE UniProtAccNum = ?
                LIMIT 1)
            """, (upac,))
        row = curs.fetchone()
        if row[0] == 1:
            result = True
        else:
            result = False
        return result

    def lookup_protein(self, upac):
        conn = self.get_conn()
        curs = conn.cursor()
        curs.execute("""
            SELECT UniProtAccNum, GeneSymbol, Name, ProcessFunction, Ions,
                Gating, InBETSE, IonChannelClassDesc, IonChannelSubClass,
                ChemblId
            FROM Protein
            WHERE UniProtAccNum = ?
            """, (upac,))
        resultset = curs.fetchone()
        if resultset == None:
            result = ['','','','','','','','','','']
        else:
            result = resultset
        return result

    def get_uniprot_accnums_by_gene_symbol(self, gene_symbol):
        conn = self.get_conn()
        curs = conn.cursor()
        curs.execute("""
            SELECT UniProtAccNum
            FROM Protein
            WHERE GeneSymbol = ?
            """, (gene_symbol,))
        resultset = curs.fetchall()
        if resultset == None:
            result = []
        else:
            result = resultset
        return resultset

    def get_gene_symbols_by_uniprot_accnum(self, upac):
        conn = self.get_conn()
        curs = conn.cursor()
        curs.execute("""
            SELECT GeneSymbol
            FROM Protein
            WHERE UniProtAccNum = ?
            """, (upac,))
        resultset = curs.fetchone()
        if resultset == None:
            result = ''
        else:
            result = resultset[0]
        return result

    def is_in_betse(self, upac):
        conn = self.get_conn()
        curs = conn.cursor()
        curs.execute("""
            SELECT InBETSE
            FROM Protein
            WHERE UniProtAccNum = ?
            """, (upac,))
        resultset = curs.fetchone()
        if resultset == None:
            result = False
        else:
            if resultset[0] == 'Y':
                result = True
            else:
                result = False
        return result

    def get_protein_uniprot_accnums(self):
        conn = self.get_conn()
        curs = conn.cursor()
        curs.execute("""
            SELECT UniProtAccNum
            FROM Protein
            """)
        resultset = curs.fetchall()
        return resultset

    def dump_protein_table(self):
        conn = self.get_conn()
        curs = conn.cursor()
        with codecs.open('protein_dump.sql', 'w', 'utf-8') as f:
            for line in conn.iterdump():
                f.write('%s\n' % line)

    def delete_proteins(self):
        conn = self.get_conn()
        curs = conn.cursor()
        curs.execute("""
            DELETE FROM Protein
            """)
        conn.commit()


    # Routines for Expression table

    def create_expression_table(self):
        conn = self.get_conn()
        curs = conn.cursor()
        curs.execute("""
            DROP TABLE IF EXISTS Expression
            """)
        curs.execute("""
            CREATE TABLE Expression(
                Id INTEGER PRIMARY KEY AUTOINCREMENT NOT NULL,
                TissueName TEXT NOT NULL,
                ProteinUniProtAccNum TEXT NOT NULL,
                ExprLevel REAL NOT NULL,
                ExprLevelQual TEXT NOT NULL,
                ExprUnits TEXT NOT NULL,
                AssayType TEXT NOT NULL,
                DatasetName TEXT NOT NULL,
                SourceDBName TEXT NOT NULL)
            """)
        conn.commit()

    def add_expression(
            self, tissue_name, upac, expr_level, expr_level_qual='',
            expr_units='', assay_type='', dataset_name='', sourcedb_name=''):
        conn = self.get_conn()
        curs = conn.cursor()
        curs.execute("""
            INSERT INTO Expression (
                TissueName, ProteinUniProtAccNum, ExprLevel, ExprLevelQual,
                ExprUnits, AssayType, DatasetName, SourceDBName)
            VALUES (?, ?, ?, ?, ?, ?, ?, ?)
            """, (
                tissue_name, upac, expr_level, expr_level_qual, expr_units,
                assay_type, dataset_name, sourcedb_name))
        conn.commit()

    def exists_expression(self, upac):
        conn = self.get_conn()
        curs = conn.cursor()
        curs.execute("""
            SELECT EXISTS(
                SELECT 1
                FROM Expression
                WHERE ProteinUniProtAccNum = ?
                LIMIT 1)
            """, (upac,))
        row = curs.fetchone()
        if row[0] == 1:
            result = True
        else:
            result = False
        return result

    def lookup_expression(self, id):
        conn = self.get_conn()
        curs = conn.cursor()
        curs.execute("""
            SELECT Id, TissueName, ProteinUniProtAccNum, ExprLevel,
                ExprLevelQual, ExprUnits, AssayType, DatasetName, SourceDBName
            FROM Expression
            WHERE Id = ?
            """, (id,))
        resultset = curs.fetchone()
        if resultset == None:
            result = ['', '', '', 0.0, '', '', '', '', '']
        else:
            result = resultset
        return result
    
    def get_expression_ids_by_uniprot_accnum(self, upac):
        conn = self.get_conn()
        curs = conn.cursor()
        curs.execute("""
            SELECT Id
            FROM Expression
            WHERE ProteinUniProtAccNum = ?
            """, (upac,))
        resultset = curs.fetchall()
        return resultset

    def get_expression_level_by_uniprot_accnum_tissue_dataset(
            self, upac, tissue, dataset):
        conn = self.get_conn()
        curs = conn.cursor()
        curs.execute("""
            SELECT ExprLevel
            FROM Expression
            WHERE ProteinUniProtAccNum = ?
                AND TissueName = ?
                AND DatasetName = ?
            """,
            (upac, tissue, dataset))
        resultset = curs.fetchall()
        return resultset

    def delete_expressions(self):
        conn = self.get_conn()
        curs = conn.cursor()
        curs.execute("""
            DELETE FROM Expression
            """)
        conn.commit()


    # Routines for DBGene table

    def create_db_gene_table(self):
        conn = self.get_conn()
        curs = conn.cursor()
        curs.execute("""
            DROP TABLE IF EXISTS DBGene
            """)
        curs.execute("""
            CREATE TABLE DBGene(
                Id INTEGER PRIMARY KEY AUTOINCREMENT NOT NULL,
                ExternalDBName TEXT NOT NULL,
                GenbankAccNum TEXT NOT NULL,
                ProbeID TEXT NOT NULL)
            """)
        conn.commit()

    def add_db_gene(self, externaldb_name, gbac, probe_id):
        conn = self.get_conn()
        curs = conn.cursor()
        curs.execute("""
            INSERT INTO DBGene(ExternalDBName, GenbankAccNum, ProbeID)
            VALUES (?, ?, ?)
            """, (externaldb_name, gbac, probe_id))
        conn.commit()


    # Routines for GenbankUniprot table

    def create_genbank_uniprot_table(self):
        conn = self.get_conn()
        curs = conn.cursor()
        curs.execute("""
            DROP TABLE IF EXISTS GenbankUniprot
            """)
        curs.execute("""
            CREATE TABLE GenbankUniprot(
                Id INTEGER PRIMARY KEY AUTOINCREMENT NOT NULL,
                GenbankAccNum TEXT NOT NULL,
                UniProtAccNum TEXT NOT NULL)
            """)
        conn.commit()

    def add_genbank_uniprot(self, gbac, upac):
        conn = self.get_conn()
        curs = conn.cursor()
        curs.execute("""
            INSERT INTO GenbankUniprot(GenbankAccNum, UniProtAccNum)
            VALUES (?, ?)
            """, (gbac, upac))
        conn.commit()

    def exists_genbank_uniprot(self, gbac, upac):
        conn = self.get_conn()
        curs = conn.cursor()
        curs.execute("""
            SELECT EXISTS(
                SELECT 1
                FROM GenbankUniprot
                WHERE GenbankAccNum=?
                    AND UniProtAccNum=?
                LIMIT 1)
            """, (gbac, upac))
        row = curs.fetchone()
        if row[0] == 1:
            result = True
        else:
            result = False
        return result

    def delete_genbank_uniprots(self):
        conn = self.get_conn()
        curs = conn.cursor()
        curs.execute("""
            DELETE FROM GenbankUniprot
            """)
        conn.commit()


    # Routines for Compound table

    def create_compound_table(self):
        conn = self.get_conn()
        curs = conn.cursor()
        curs.execute("""
            DROP TABLE IF EXISTS Compound
            """)
        curs.execute("""
            CREATE TABLE Compound(
                Id INTEGER PRIMARY KEY AUTOINCREMENT NOT NULL,
                SMILES TEXT NOT NULL,
                InChI TEXT NOT NULL,
                Name TEXT NOT NULL,
                ChemblId TEXT NOT NULL UNIQUE,
                Synonyms TEXT NOT NULL,
                ApprovalStatus TEXT NOT NULL,
                FirstApprovalYear TEXT NOT NULL,
                SourceDBName NOT NULL)
            """)
        conn.commit()

    def add_compound(
            self, name, smiles='', inchi='', chembl_id='', synonyms='',
            approval_status='', first_approval_year='', sourcedb_name=''):
        conn = self.get_conn()
        curs = conn.cursor()
        curs.execute("""
            INSERT INTO Compound(
                SMILES, InChI, Name, ChemblId, Synonyms, ApprovalStatus,
                FirstApprovalYear, SourceDBName)
            VALUES (?, ?, ?, ?, ?, ?, ?, ?)
            """, (
                smiles, inchi, name, chembl_id, synonyms, approval_status,
                first_approval_year, sourcedb_name))
        conn.commit()

    def exists_compound_by_chembl_id(self, chembl_id):
        conn = self.get_conn()
        curs = conn.cursor()
        expr_list = []
        curs.execute("""
            SELECT EXISTS(
                SELECT 1
                FROM Compound
                WHERE ChemblId = ?
                LIMIT 1)
            """, (chembl_id,))
        row = curs.fetchone()
        if row[0] == 1:
            result = True
        else:
            result = False
        return result

    def lookup_compound(self, id):
        conn = self.get_conn()
        curs = conn.cursor()
        curs.execute("""
            SELECT Id, SMILES, InChI, Name, ChemblId, Synonyms, ApprovalStatus,
                FirstApprovalYear, SourceDBName
            FROM Compound
            WHERE Id = ?
            """, (id,))
        resultset = curs.fetchone()
        if resultset == None:
            result = ['', '', '', '', '', '', '', '', '']
        else:
            result = resultset
        return result

    def get_compound_id_by_chembl_id(self, chembl_id):
        conn = self.get_conn()
        curs = conn.cursor()
        curs.execute("""
            SELECT Id
            FROM Compound
            WHERE ChemblId = ?
            """, (chembl_id,))
        resultset = curs.fetchone()
        if resultset == None:
            result = ''
        else:
            result = resultset[0]
        return result

    def get_compound_id_by_compound_name(self, compound_name):
        conn = self.get_conn()
        curs = conn.cursor()
        curs.execute("""
            SELECT Id
            FROM Compound
            WHERE Name = ?
            """, (compound_name,))
        resultset = curs.fetchone()
        if resultset == None:
            result = ''
        else:
            result = resultset[0]
        return result

    def delete_compounds(self):
        conn = self.get_conn()
        curs = conn.cursor()
        curs.execute("""
            DELETE FROM Compound
            """)
        conn.commit()


    # Routines for Interaction table

    def create_interaction_table(self):
        conn = self.get_conn()
        curs = conn.cursor()
        curs.execute("""
            DROP TABLE IF EXISTS Interaction
            """)
        curs.execute("""
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
                SourceDBName TEXT NOT NULL)
            """)
        conn.commit()

    def add_interaction(
            self, target_upac, compound_id, action_type='', action_desc='',
            strength=0, strength_units='', assay_type='', chembl_id='',
            sourcedb_name=''):
        conn = self.get_conn()
        curs = conn.cursor()
        curs.execute("""
            INSERT INTO Interaction(
                TargetUniProtAccNum, CompoundId, ActionType, ActionDesc,
                Strength, StrengthUnits, AssayType, ChemblId, SourceDBName)
            VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?)
            """, (
                target_upac, compound_id, action_type, action_desc, strength,
                strength_units, assay_type, chembl_id, sourcedb_name))
        conn.commit()

    def exists_interaction_uniprot_accnum_compound_id(self, upac, compound_id):
        conn = self.get_conn()
        curs = conn.cursor()
        curs.execute("""
            SELECT EXISTS(
                SELECT 1
                FROM Interaction
                WHERE TargetUniProtAccNum=?
                    AND CompoundId=?
                LIMIT 1)
            """, (upac, compound_id))
        row = curs.fetchone()
        if row[0] == 1:
            result = True
        else:
            result = False
        return result

    def exists_interaction_by_uniprot_accnum(self, upac):
        conn = self.get_conn()
        curs = conn.cursor()
        curs.execute("""
            SELECT EXISTS(
                SELECT 1
                FROM Interaction
                WHERE TargetUniProtAccNum = ?
                LIMIT 1)
            """, (upac,))
        row = curs.fetchone()
        if row[0] == 1:
            result = True
        else:
            result = False
        return result

    def get_interaction_ids_by_uniprot_accnum(self, upac):
        conn = self.get_conn()
        curs = conn.cursor()
        curs.execute("""
            SELECT Id
            FROM Interaction
            WHERE TargetUniProtAccNum = ?
            """, (upac,))
        resultset = curs.fetchall()
        return resultset

    def lookup_interaction(self, id):
        conn = self.get_conn()
        curs = conn.cursor()
        curs.execute("""
            SELECT Id, TargetUniProtAccNum, CompoundId, ActionType, ActionDesc,
                Strength, StrengthUnits, AssayType, ChemblId, SourceDBName
            FROM Interaction
            WHERE Id = ?
            """, (id,))
        resultset = curs.fetchone()
        if resultset == None:
            result = ['', '', '', '', '', '', '', '', '', '']
        else:
            result = resultset[0]
        return resultset

    def delete_interactions(self):
        conn = self.get_conn()
        curs = conn.cursor()
        curs.execute("""
            DELETE FROM Interaction
            """)
        conn.commit()


    # Routines for DBTissue table

    def create_db_tissue_table(self):
        conn = self.get_conn()
        curs = conn.cursor()
        curs.execute("""
            DROP TABLE IF EXISTS DBTissue
            """)
        curs.execute("""
            CREATE TABLE DBTissue(
                Id INTEGER PRIMARY KEY AUTOINCREMENT NOT NULL,
                ExternalDBName TEXT NOT NULL,
                TissueName TEXT NOT NULL,
                DBEquivalentTissueName TEXT NOT NULL)
            """)
        conn.commit()

    def add_db_tissue(
            self, externaldb_name, tissue_name, db_equivalent_tissue_name):
        conn = self.get_conn()
        curs = conn.cursor()
        curs.execute("""
            INSERT INTO DBTissue(
                ExternalDBName, TissueName, DBEquivalentTissueName)
            VALUES (?, ?, ?)
            """, (externaldb_name, tissue_name, db_equivalent_tissue_name))
        conn.commit()
    
    def get_tissue_name_by_db_tissue_name(
            self, externaldb_name, db_equivalent_tissue_name):
        conn = self.get_conn()
        curs = conn.cursor()
        curs.execute("""
            SELECT TissueName
            FROM DBTissue
            WHERE ExternalDBName=?
                AND DBEquivalentTissueName=?
            """, (externaldb_name, db_equivalent_tissue_name))
        resultset = curs.fetchone()
        if resultset == None:
            result = ''
        else:
            result = resultset[0]
        return result

    def delete_db_tissues(self):
        conn = self.get_conn()
        curs = conn.cursor()
        curs.execute("""
            DELETE FROM DBTissue
            """)
        conn.commit()


    # Routines for Specificity table

    def create_specificity_table(self):
        conn = self.get_conn()
        curs = conn.cursor()
        curs.execute("""
            DROP TABLE IF EXISTS Specificity
            """)
        curs.execute("""
            CREATE TABLE Specificity(
                Id INTEGER PRIMARY KEY AUTOINCREMENT NOT NULL,
                TissueName TEXT NOT NULL,
                UniProtAccNum TEXT NOT NULL,
                SpecificityScore TEXT NOT NULL)
            """)
        conn.commit()

    def add_specificity(self, tissue_name, upac, specificity_score):
        conn = self.get_conn()
        curs = conn.cursor()
        curs.execute("""
            INSERT INTO Specificity(
                TissueName, UniProtAccNum, SpecificityScore)
            VALUES (?, ?, ?)
            """, (tissue_name, upac, specificity_score))
        conn.commit()

    def delete_specificities(self):
        conn = self.get_conn()
        curs = conn.cursor()
        curs.execute("""
            DELETE FROM Specificity
            """)
        conn.commit()


    # Routines for GoTerm table

    def create_go_term_table(self):
        conn = self.get_conn()
        curs = conn.cursor()
        curs.execute("""
            DROP TABLE IF EXISTS GoTerm
            """)
        curs.execute("""
            CREATE TABLE GoTerm(
                Id INTEGER PRIMARY KEY AUTOINCREMENT NOT NULL,
                UniProtAccNum TEXT NOT NULL,
                GoId TEXT NOT NULL)
            """)
        conn.commit()

    def add_go_term(self, upac, goid):
        conn = self.get_conn()
        curs = conn.cursor()
        curs.execute("""
            INSERT INTO GoTerm(UniProtAccNum, GoId)
            VALUES (?, ?)
            """, (upac, goid))
        conn.commit()

    def delete_go_terms(self):
        conn = self.get_conn()
        curs = conn.cursor()
        curs.execute("""
            DELETE FROM GoTerm
            """)
        conn.commit()


    # Routines for PDB table
    
    def create_pdb_table(self):
        conn = self.get_conn()
        curs = conn.cursor()
        curs.execute("""
            DROP TABLE IF EXISTS PDB
            """)
        curs.execute("""
            CREATE TABLE PDB(
                Id INTEGER PRIMARY KEY AUTOINCREMENT NOT NULL,
                UniProtAccNum TEXT NOT NULL,
                PDBID TEXT NOT NULL,
                HumanOnly TEXT NOT NULL)
            """)
        conn.commit()


    # Other database routines

    def vacuum_db(self):
        conn = self.get_conn()
        curs = conn.cursor()
        curs.execute("""
            VACUUM
            """)
        conn.commit()

    def dump_db(self):
        conn = self.get_conn()
        with codecs.open('icdb_dump.sql', 'w', 'utf-8') as f:
            for line in conn.iterdump():
                f.write('%s\n' % line)

    def print_db_stats(self):
        conn = self.get_conn()
        curs = conn.cursor()
        curs.execute("""
            SELECT *
            FROM sqlite_master
            WHERE type = 'table'
            """)
        table_name_list = []
        for table_info in curs:
            table_name = table_info[1]
            table_name_list.append(table_name)
        for table_name in table_name_list:
            curs.execute("""
                SELECT COUNT(*)
                FROM %s
                """ % (table_name,))
            row_count = curs.fetchone()[0]
            print table_name, row_count
        useful_data_count = 0
        upac_record_list = self.get_protein_uniprot_accnums()
        for upac_record in upac_record_list:
            upac = upac_record[0]
            interaction_id_list = \
                self.get_interaction_ids_by_uniprot_accnum(upac)
            expression_id_list = \
                self.get_expression_ids_by_uniprot_accnum(upac)
            if len(interaction_id_list) > 0 and len(expression_id_list) > 0:
                useful_data_count += 1
        print 'Protein with interaction and expression %d' % useful_data_count
        upac_record_list = self.get_protein_uniprot_accnums()
        in_betse_count = 0
        in_betse_with_expr_count = 0
        list_in_betse_with_expr = []
        list_in_betse_no_expr = []
        for upac_record in upac_record_list:
            upac = upac_record[0]
            protein_record = self.lookup_protein(upac)
            in_betse = protein_record[6]
            if in_betse == 'Y':
                in_betse_count += 1
                if self.exists_expression(upac):
                    in_betse_with_expr_count += 1
                    list_in_betse_with_expr.append(
                        self.get_gene_symbols_by_uniprot_accnum(upac))
                else:
                    list_in_betse_no_expr.append(
                        self.get_gene_symbols_by_uniprot_accnum(upac))
        print 'Number of proteins in betse = %d' % in_betse_count
        print '...also with with expression = %d' % in_betse_with_expr_count
