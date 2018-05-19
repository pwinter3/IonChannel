import codecs
import sqlite3


class IonChannelDatabase(object):

    def __init__(self, path_db):
        self._path_db = path_db
        self._conn = sqlite3.connect(path_db)
        self._conn.text_factory = str

    def _get_conn(self):
        return self._conn
    
    
    # Routines for executing multiple queries in one transaction
    # Do not work

    def begin_transaction(self):
        conn = self._get_conn()
        conn.execute('BEGIN TRANSACTION')

    def end_transaction(self):
        conn = self._get_conn()
        conn.execute('COMMIT')


    # Routines for ChannelSuperClass table

    def create_channel_super_class_table(self):
        conn = self._get_conn()
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
        conn = self._get_conn()
        curs = conn.cursor()
        curs.execute("""
            INSERT INTO ChannelSuperClass(Name)
            VALUES (?)
            """, (name,))
        conn.commit()

    def exists_channel_super_class(self, name):
        conn = self._get_conn()
        curs = conn.cursor()
        curs.execute("""
            SELECT EXISTS(
                SELECT 1
                FROM ChannelSuperClass
                WHERE Name = ?
                LIMIT 1)
            """, (name,))
        row = curs.fetchone()
        return bool(row[0] == 1)

    def delete_channel_super_classes(self):
        conn = self._get_conn()
        curs = conn.cursor()
        curs.execute("""
            DELETE FROM ChannelSuperClass
            """)
        conn.commit()


    # Routines for ChannelClass table

    def create_channel_class_table(self):
        conn = self._get_conn()
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
        conn = self._get_conn()
        curs = conn.cursor()
        curs.execute("""
            INSERT INTO ChannelClass(Name, SuperClass)
            VALUES (?, ?)
            """, (name, super_class))
        conn.commit()

    def exists_channel_class(self, name):
        conn = self._get_conn()
        curs = conn.cursor()
        curs.execute("""
            SELECT EXISTS(
                SELECT 1
                FROM ChannelClass
                WHERE Name = ?
                LIMIT 1)
            """, (name,))
        row = curs.fetchone()
        return bool(row[0] == 1)

    def delete_channel_classes(self):
        conn = self._get_conn()
        curs = conn.cursor()
        curs.execute("""
            DELETE FROM ChannelClass
            """)
        conn.commit()


    # Routines for ChannelSubClass table

    def create_channel_sub_class_table(self):
        conn = self._get_conn()
        curs = conn.cursor()
        curs.execute("""
            DROP TABLE IF EXISTS ChannelSubClass
            """)
        curs.execute("""
            CREATE TABLE ChannelSubClass(
                Name TEXT PRIMARY KEY NOT NULL,
                Class TEXT NOT NULL,
                ChannelpediaText TEXT NOT NULL,
                ChannelpediaURL TEXT NOT NULL,
                Subfamily TEXT NOT NULL)
            """)
        conn.commit()

    def add_channel_sub_class(
            self, name, channel_class='', channelpedia_text='',
            channelpedia_url='', subfamily=''):
        conn = self._get_conn()
        conn.text_factory = str
        curs = conn.cursor()
        curs.execute("""
            INSERT INTO ChannelSubClass(
                Name, Class, ChannelpediaText, ChannelpediaURL, Subfamily)
            VALUES (?, ?, ?, ?, ?)
            """, (
                name, channel_class, channelpedia_text, channelpedia_url,
                subfamily))
        conn.commit()

    def exists_channel_sub_class(self, name):
        conn = self._get_conn()
        curs = conn.cursor()
        curs.execute("""
            SELECT EXISTS(
                SELECT 1
                FROM ChannelSubClass
                WHERE Name = ?
                LIMIT 1)
            """, (name,))
        row = curs.fetchone()
        return bool(row[0] == 1)

    def get_channelpedia_info(self, name):
        conn = self._get_conn()
        conn.row_factory = sqlite3.Row
        curs = conn.cursor()
        curs.execute("""
            SELECT ChannelpediaText, ChannelpediaURL
            FROM ChannelSubClass
            WHERE Name = ?
            """, (name,))
        row = curs.fetchone()
        return row

    def delete_channel_sub_classes(self):
        conn = self._get_conn()
        curs = conn.cursor()
        curs.execute("""
            DELETE FROM ChannelSubClass
            """)
        conn.commit()


    # Routines for ExternalDB table

    def create_external_db_table(self):
        conn = self._get_conn()
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
        conn = self._get_conn()
        curs = conn.cursor()
        curs.execute("""
            INSERT INTO ExternalDB(Name, URL)
            VALUES (?, ?)
            """, (name, url))
        conn.commit()

    def delete_external_dbs(self):
        conn = self._get_conn()
        curs = conn.cursor()
        curs.execute("""
            DELETE FROM ExternalDB
            """)
        conn.commit()


    # Routines for Tissue table

    def create_tissue_table(self):
        conn = self._get_conn()
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
        conn = self._get_conn()
        curs = conn.cursor()
        curs.execute("""
            INSERT INTO Tissue(Name, BTOId, ParentTissueName)
            VALUES (?, ?, ?)
            """, (name, bto_id, parent_tissue_name))
        conn.commit()

    def get_tissue_names(self):
        conn = self._get_conn()
        curs = conn.cursor()
        curs.execute("""
            SELECT Name
            FROM Tissue
            """)
        row_list = curs.fetchall()
        return [elem[0] for elem in row_list]

    def delete_tissues(self):
        conn = self._get_conn()
        curs = conn.cursor()
        curs.execute("""
            DELETE FROM Tissue
            """)
        conn.commit()


    # Routines for Protein table

    def create_protein_table(self):
        conn = self._get_conn()
        curs = conn.cursor()
        curs.execute("""
            DROP TABLE IF EXISTS Protein
            """)
        curs.execute("""
            CREATE TABLE Protein(
                UniProtAccNum TEXT PRIMARY KEY NOT NULL,
                GeneSymbol TEXT NOT NULL,
                Name TEXT NOT NULL,
                Ions TEXT NOT NULL,
                Gating TEXT NOT NULL,
                InBETSE TEXT NOT NULL,
                IonChannelSubClass TEXT NOT NULL,
                ChemblId TEXT NOT NULL)
            """)
        conn.commit()

    def add_protein(
            self, upac, gene_symbol='', name='', ions='', gating='',
            in_betse='N', ion_channel_sub_class='', chembl_id=''):
        conn = self._get_conn()
        curs = conn.cursor()
        curs.execute("""
            INSERT INTO Protein(
                UniProtAccNum, GeneSymbol, Name, Ions, Gating, InBETSE,
                IonChannelSubClass, ChemblId)
            VALUES (?, ?, ?, ?, ?, ?, ?, ?)
            """, (
                upac, gene_symbol, name, ions, gating, in_betse,
                ion_channel_sub_class, chembl_id))
        conn.commit()

    def update_protein_sub_class(self, upac, ion_channel_sub_class):
        conn = self._get_conn()
        curs = conn.cursor()
        curs.execute("""
            UPDATE Protein
            SET IonChannelSubClass = ?
            WHERE UniProtAccNum = ?
            """, (ion_channel_sub_class, upac))
        conn.commit()

    def update_protein_name(self, upac, name):
        conn = self._get_conn()
        curs = conn.cursor()
        curs.execute("""
            UPDATE Protein
            SET Name = ?
            WHERE UniProtAccNum = ?
            """, (name, upac))
        conn.commit()

    def update_protein_chembl_id(self, upac, chembl_id):
        conn = self._get_conn()
        curs = conn.cursor()
        curs.execute("""
            UPDATE Protein
            SET ChemblId = ?
            WHERE UniProtAccNum = ?
            """, (chembl_id, upac))
        conn.commit()

    def update_protein_gene_symbol(self, upac, gene_symbol):
        conn = self._get_conn()
        curs = conn.cursor()
        curs.execute("""
            UPDATE Protein
            SET GeneSymbol = ?
            WHERE UniProtAccNum = ?
            """, (gene_symbol, upac))
        conn.commit()

    def update_protein_ions(self, upac, ions):
        conn = self._get_conn()
        curs = conn.cursor()
        curs.execute("""
            UPDATE Protein
            SET Ions = ?
            WHERE UniProtAccNum = ?
            """, (ions, upac))
        conn.commit()

    def update_protein_gating(self, upac, gating):
        conn = self._get_conn()
        curs = conn.cursor()
        curs.execute("""
            UPDATE Protein SET Gating = ?
            WHERE UniProtAccNum = ?
            """, (gating, upac))
        conn.commit()

    def update_protein_in_betse(self, upac, in_betse):
        conn = self._get_conn()
        curs = conn.cursor()
        curs.execute("""
            UPDATE Protein
            SET InBETSE = ?
            WHERE UniProtAccNum = ?
            """, (in_betse, upac))
        conn.commit()

    def update_protein_ion_channel_sub_class(
            self, upac, ion_channel_sub_class):
        conn = self._get_conn()
        curs = conn.cursor()
        curs.execute("""
            UPDATE Protein
            SET IonChannelSubClass = ?
            WHERE UniProtAccNum = ?
            """, (ion_channel_sub_class, upac))
        conn.commit()

    def exists_protein(self, upac):
        conn = self._get_conn()
        curs = conn.cursor()
        curs.execute("""
            SELECT EXISTS(
                SELECT 1
                FROM Protein
                WHERE UniProtAccNum = ?
                LIMIT 1)
            """, (upac,))
        row = curs.fetchone()
        return bool(row[0] == 1)

    def lookup_protein(self, upac):
        conn = self._get_conn()
        conn.row_factory = sqlite3.Row
        curs = conn.cursor()
        curs.execute("""
            SELECT UniProtAccNum, GeneSymbol, Name, Ions, Gating, InBETSE,
                IonChannelSubClass, ChemblId
            FROM Protein
            WHERE UniProtAccNum = ?
            """, (upac,))
        row = curs.fetchone()
        return row

    def get_uniprot_accnums_by_gene_symbol(self, gene_symbol):
        conn = self._get_conn()
        curs = conn.cursor()
        curs.execute("""
            SELECT UniProtAccNum
            FROM Protein
            WHERE GeneSymbol = ?
            """, (gene_symbol,))
        row_list = curs.fetchall()
        return [elem[0] for elem in row_list]

    def get_gene_symbols_by_uniprot_accnum(self, upac):
        conn = self._get_conn()
        curs = conn.cursor()
        curs.execute("""
            SELECT GeneSymbol
            FROM Protein
            WHERE UniProtAccNum = ?
            """, (upac,))
        row_list = curs.fetchall()
        return [elem[0] for elem in row_list]

    def is_in_betse(self, upac):
        conn = self._get_conn()
        curs = conn.cursor()
        curs.execute("""
            SELECT InBETSE
            FROM Protein
            WHERE UniProtAccNum = ?
            """, (upac,))
        row = curs.fetchone()
        return bool(row[0] == 'Y')

    def get_protein_uniprot_accnums(self):
        conn = self._get_conn()
        curs = conn.cursor()
        curs.execute("""
            SELECT UniProtAccNum
            FROM Protein
            """)
        row_list = curs.fetchall()
        return [elem[0] for elem in row_list]

    def get_in_betse_or_expr_threshold_protein(
            self, tissue_list, threshold, include_betse=True,
            specificity_option='comprehensive'):
        conn = self._get_conn()
        curs = conn.cursor()
        if include_betse:
            sql = """
                SELECT Protein.UniProtAccNum, Protein.Name, GeneSymbol
                FROM Protein
                INNER JOIN Expression
                    ON Protein.UniProtAccNum = Expression.ProteinUniProtAccNum
                INNER JOIN GoTerm
                    ON Protein.UniProtAccNum = GoTerm.UniProtAccNum
                WHERE (
                    Protein.InBETSE = "Y" OR (
                        GoTerm.GoId IN (
                                "GO:0001518","GO:0005891","GO:0005892",
                                "GO:0008076","GO:0008328","GO:0016935",
                                "GO:0017071","GO:0017146","GO:0032281",
                                "GO:0032983","GO:0034702","GO:0034703",
                                "GO:0034704","GO:0034705","GO:0034706",
                                "GO:0034707","GO:0036128","GO:0071193",
                                "GO:0098855","GO:1902937","GO:1904602",
                                "GO:1990246","GO:1990425","GO:1990454")
                            AND Expression.ExprLevel > ?
                            AND Expression.TissueName IN (
                """ + ','.join(['?'] * len(tissue_list)) + ')))'
        else:
            sql = """
                SELECT Protein.UniProtAccNum, Protein.Name, GeneSymbol
                FROM Protein
                INNER JOIN Expression
                    ON Protein.UniProtAccNum = Expression.ProteinUniProtAccNum
                INNER JOIN GoTerm
                    ON Protein.UniProtAccNum = GoTerm.UniProtAccNum
                WHERE (
                    Expression.ExprLevel > ?
                        AND GoTerm.GoId IN (
                            "GO:0001518","GO:0005891","GO:0005892",
                            "GO:0008076","GO:0008328","GO:0016935",
                            "GO:0017071","GO:0017146","GO:0032281",
                            "GO:0032983","GO:0034702","GO:0034703",
                            "GO:0034704","GO:0034705","GO:0034706",
                            "GO:0034707","GO:0036128","GO:0071193",
                            "GO:0098855","GO:1902937","GO:1904602",
                            "GO:1990246","GO:1990425","GO:1990454")
                        AND Expression.TissueName IN (
                """ + ','.join(['?'] * len(tissue_list)) + '))'
        curs.execute(sql, (threshold,) + tuple(tissue_list))
        resultset = curs.fetchall()
        resultset = list(set(resultset)) # Remove redundant records
        return resultset

    def delete_proteins(self):
        conn = self._get_conn()
        curs = conn.cursor()
        curs.execute("""
            DELETE FROM Protein
            """)
        conn.commit()


    # Routines for Expression table

    def create_expression_table(self):
        conn = self._get_conn()
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
            self, tissue_name, upac, expr_level=0., expr_level_qual='',
            expr_units='', assay_type='', dataset_name='', source_db_name='',
            commit=True):
        conn = self._get_conn()
        curs = conn.cursor()
        curs.execute("""
            INSERT INTO Expression (
                TissueName, ProteinUniProtAccNum, ExprLevel, ExprLevelQual,
                ExprUnits, AssayType, DatasetName, SourceDBName)
            VALUES (?, ?, ?, ?, ?, ?, ?, ?)
            """, (
                tissue_name, upac, expr_level, expr_level_qual, expr_units,
                assay_type, dataset_name, source_db_name))
        if commit:
            conn.commit()

    def exists_expression(self, upac):
        conn = self._get_conn()
        curs = conn.cursor()
        curs.execute("""
            SELECT EXISTS(
                SELECT 1
                FROM Expression
                WHERE ProteinUniProtAccNum = ?
                LIMIT 1)
            """, (upac,))
        row = curs.fetchone()
        return bool(row[0] == 1)

    def lookup_expression(self, expression_id):
        conn = self._get_conn()
        conn.row_factory = sqlite3.Row
        curs = conn.cursor()
        curs.execute("""
            SELECT Id, TissueName, ProteinUniProtAccNum, ExprLevel,
                ExprLevelQual, ExprUnits, AssayType, DatasetName, SourceDBName
            FROM Expression
            WHERE Id = ?
            """, (expression_id,))
        row = curs.fetchone()
        return row

    def get_expression_ids_by_uniprot_accnum(self, upac):
        conn = self._get_conn()
        curs = conn.cursor()
        curs.execute("""
            SELECT Id
            FROM Expression
            WHERE ProteinUniProtAccNum = ?
            """, (upac,))
        row_list = curs.fetchall()
        return [elem[0] for elem in row_list]

    def get_expression_level_by_uniprot_accnum_tissue_dataset(
            self, upac, tissue_name, dataset_name):
        conn = self._get_conn()
        curs = conn.cursor()
        curs.execute("""
            SELECT ExprLevel
            FROM Expression
            WHERE ProteinUniProtAccNum = ?
                AND TissueName = ?
                AND DatasetName = ?
            """, (upac, tissue_name, dataset_name))
        row_list = curs.fetchall()
        return [elem[0] for elem in row_list]

    def get_expression_level_qual_by_uniprot_accnum_tissue_dataset(
            self, upac, tissue_name, dataset_name):
        conn = self._get_conn()
        curs = conn.cursor()
        curs.execute("""
            SELECT ExprLevelQual
            FROM Expression
            WHERE ProteinUniProtAccNum = ?
                AND TissueName = ?
                AND DatasetName = ?
            """, (upac, tissue_name, dataset_name))
        row_list = curs.fetchall()
        return [elem[0] for elem in row_list]

    def delete_expressions(self):
        conn = self._get_conn()
        curs = conn.cursor()
        curs.execute("""
            DELETE FROM Expression
            """)
        conn.commit()


    # Routines for DBGene table

    def create_db_gene_table(self):
        conn = self._get_conn()
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

    def add_db_gene(self, external_db_name, gbac, probe_id):
        conn = self._get_conn()
        curs = conn.cursor()
        curs.execute("""
            INSERT INTO DBGene(ExternalDBName, GenbankAccNum, ProbeID)
            VALUES (?, ?, ?)
            """, (external_db_name, gbac, probe_id))
        conn.commit()


    # Routines for GenbankUniprot table
    # Table is not being used, should remove this code

    def create_genbank_uniprot_table(self):
        conn = self._get_conn()
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
        conn = self._get_conn()
        curs = conn.cursor()
        curs.execute("""
            INSERT INTO GenbankUniprot(GenbankAccNum, UniProtAccNum)
            VALUES (?, ?)
            """, (gbac, upac))
        conn.commit()

    def exists_genbank_uniprot(self, gbac, upac):
        conn = self._get_conn()
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
        return bool(row[0] == 1)

    def delete_genbank_uniprots(self):
        conn = self._get_conn()
        curs = conn.cursor()
        curs.execute("""
            DELETE FROM GenbankUniprot
            """)
        conn.commit()


    # Routines for Compound table

    def create_compound_table(self):
        conn = self._get_conn()
        curs = conn.cursor()
        curs.execute("""
            DROP TABLE IF EXISTS Compound
            """)
        curs.execute("""
            CREATE TABLE Compound(
                Id INTEGER PRIMARY KEY AUTOINCREMENT NOT NULL,
                CompoundType TEXT NOT NULL,
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
            self, name, compound_type='', smiles='', inchi='', chembl_id='',
            synonyms='', approval_status='', first_approval_year='',
            source_db_name=''):
        conn = self._get_conn()
        curs = conn.cursor()
        curs.execute("""
            INSERT INTO Compound(
                CompoundType, SMILES, InChI, Name, ChemblId, Synonyms,
                ApprovalStatus, FirstApprovalYear, SourceDBName)
            VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?)
            """, (
                compound_type, smiles, inchi, name, chembl_id, synonyms,
                approval_status, first_approval_year, source_db_name))
        conn.commit()

    def exists_compound_by_chembl_id(self, chembl_id):
        conn = self._get_conn()
        curs = conn.cursor()
        curs.execute("""
            SELECT EXISTS(
                SELECT 1
                FROM Compound
                WHERE ChemblId = ?
                LIMIT 1)
            """, (chembl_id,))
        row = curs.fetchone()
        return bool(row[0] == 1)

    def lookup_compound(self, compound_id):
        conn = self._get_conn()
        conn.row_factory = sqlite3.Row
        curs = conn.cursor()
        curs.execute("""
            SELECT Id, CompoundType, SMILES, InChI, Name, ChemblId, Synonyms,
                ApprovalStatus, FirstApprovalYear, SourceDBName
            FROM Compound
            WHERE Id = ?
            """, (compound_id,))
        row = curs.fetchone()
        return row

    def get_compound_id_by_chembl_id(self, chembl_id):
        conn = self._get_conn()
        curs = conn.cursor()
        curs.execute("""
            SELECT Id
            FROM Compound
            WHERE ChemblId = ?
            """, (chembl_id,))
        row = curs.fetchone()
        return row[0]

    def get_compound_ids_by_compound_name(self, compound_name):
        conn = self._get_conn()
        curs = conn.cursor()
        curs.execute("""
            SELECT Id
            FROM Compound
            WHERE Name = ?
            """, (compound_name,))
        row_list = curs.fetchall()
        return [elem[0] for elem in row_list]

    def delete_compounds(self):
        conn = self._get_conn()
        curs = conn.cursor()
        curs.execute("""
            DELETE FROM Compound
            """)
        conn.commit()


    # Routines for Interaction table

    def create_interaction_table(self):
        conn = self._get_conn()
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
                AssayValue REAL NOT NULL,
                AssayUnits TEXT NOT NULL,
                AssayStandardType TEXT NOT NULL,
                AssayType TEXT NOT NULL,
                AssayChemblId TEXT NOT NULL,
                SourceDBName TEXT NOT NULL)
            """)
        conn.commit()

    def add_interaction(
            self, target_upac, compound_id, action_type='', assay_value=0,
            assay_units='', assay_standard_type='', assay_type='', 
            assay_chembl_id='', source_db_name=''):
        conn = self._get_conn()
        curs = conn.cursor()
        curs.execute("""
            INSERT INTO Interaction(
                TargetUniProtAccNum, CompoundId, ActionType, AssayValue,
                AssayUnits, AssayStandardType, AssayType, AssayChemblId,
                SourceDBName)
            VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?)
            """, (
                target_upac, compound_id, action_type, assay_value,
                assay_units, assay_standard_type, assay_type, assay_chembl_id,
                source_db_name))
        conn.commit()

    def exists_interaction_uniprot_accnum_compound_id(self, upac, compound_id):
        conn = self._get_conn()
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
        return bool(row[0] == 1)

    def exists_interaction_by_uniprot_accnum(self, upac):
        conn = self._get_conn()
        curs = conn.cursor()
        curs.execute("""
            SELECT EXISTS(
                SELECT 1
                FROM Interaction
                WHERE TargetUniProtAccNum = ?
                LIMIT 1)
            """, (upac,))
        row = curs.fetchone()
        return bool(row[0] == 1)

    def get_interaction_ids_by_uniprot_accnum(
            self, upac, type=None, unit=None):
        conn = self._get_conn()
        curs = conn.cursor()
        if type != None and unit != None:
            curs.execute('''
                SELECT Id
                FROM Interaction
                WHERE TargetUniProtAccNum = ?
                    AND AssayStandardType = ?
                    AND AssayUnits = ?
                ''', (upac, type, unit))
        elif type != None:
            curs.execute('''
                SELECT Id
                FROM Interaction
                WHERE TargetUniProtAccNum = ?
                    AND AssayStandardType = ?
                ''', (upac, type))
        elif unit != None:
            curs.execute('''
                SELECT Id
                FROM Interaction
                WHERE TargetUniProtAccNum = ?
                    AND AssayUnits = ?
                ''', (upac, unit))
        else:
            curs.execute('''
                SELECT Id
                FROM Interaction
                WHERE TargetUniProtAccNum = ?
                ''', (upac,))
        row_list = curs.fetchall()
        return [elem[0] for elem in row_list]

    def lookup_interaction(self, interaction_id):
        conn = self._get_conn()
        conn.row_factory = sqlite3.Row
        curs = conn.cursor()
        curs.execute("""
            SELECT Id, TargetUniProtAccNum, CompoundId, ActionType,
                AssayValue, AssayUnits, AssayStandardType, AssayChemblId,
                SourceDBName
            FROM Interaction
            WHERE Id = ?
            """, (interaction_id,))
        row = curs.fetchone()
        return row

    def delete_interactions(self):
        conn = self._get_conn()
        curs = conn.cursor()
        curs.execute("""
            DELETE FROM Interaction
            """)
        conn.commit()


    # Routines for DBTissue table

    def create_db_tissue_table(self):
        conn = self._get_conn()
        curs = conn.cursor()
        curs.execute("""
            DROP TABLE IF EXISTS DBTissue
            """)
        curs.execute("""
            CREATE TABLE DBTissue(
                Id INTEGER PRIMARY KEY AUTOINCREMENT NOT NULL,
                ExternalDBName TEXT NOT NULL,
                TissueName TEXT NOT NULL,
                DBEquivalentTissueName TEXT NOT NULL,
                CONSTRAINT UniCon UNIQUE (
                    ExternalDBName, DBEquivalentTissueName))
            """)
        conn.commit()

    def add_db_tissue(
            self, external_db_name, tissue_name, db_equivalent_tissue_name):
        conn = self._get_conn()
        curs = conn.cursor()
        curs.execute("""
            INSERT INTO DBTissue(
                ExternalDBName, TissueName, DBEquivalentTissueName)
            VALUES (?, ?, ?)
            """, (external_db_name, tissue_name, db_equivalent_tissue_name))
        conn.commit()

    def get_tissue_name_by_db_tissue_name(
            self, external_db_name, db_equivalent_tissue_name):
        conn = self._get_conn()
        curs = conn.cursor()
        curs.execute("""
            SELECT TissueName
            FROM DBTissue
            WHERE ExternalDBName = ?
                AND DBEquivalentTissueName = ?
            """, (external_db_name, db_equivalent_tissue_name))
        row = curs.fetchone()
        if row is not None:
            tissue_name = row[0]
        else:
            tissue_name = ''
        return tissue_name

    def delete_db_tissues(self):
        conn = self._get_conn()
        curs = conn.cursor()
        curs.execute("""
            DELETE FROM DBTissue
            """)
        conn.commit()


    # Routines for Specificity table

    def create_specificity_table(self):
        conn = self._get_conn()
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
        conn = self._get_conn()
        curs = conn.cursor()
        curs.execute("""
            INSERT INTO Specificity(
                TissueName, UniProtAccNum, SpecificityScore)
            VALUES (?, ?, ?)
            """, (tissue_name, upac, specificity_score))
        conn.commit()

    def get_specificity(self, tissue_name, upac):
        conn = self._get_conn()
        conn.row_factory = sqlite3.Row
        curs = conn.cursor()
        curs.execute("""
            SELECT SpecificityScore
            FROM Specificity
            WHERE TissueName = ?
                AND UniProtAccNum = ?
            """, (tissue_name, upac))
        row = curs.fetchone()
        value = row[0]
        return value

    def delete_specificities(self):
        conn = self._get_conn()
        curs = conn.cursor()
        curs.execute("""
            DELETE FROM Specificity
            """)
        conn.commit()


    # Routines for GoTerm table

    def create_go_term_table(self):
        conn = self._get_conn()
        curs = conn.cursor()
        curs.execute("""
            DROP TABLE IF EXISTS GoTerm
            """)
        curs.execute("""
            CREATE TABLE GoTerm(
                Id INTEGER PRIMARY KEY AUTOINCREMENT NOT NULL,
                UniProtAccNum TEXT NOT NULL,
                GoId TEXT NOT NULL,
                Name TEXT NOT NULL,
                Qualifier TEXT NOT NULL,
                Aspect TEXT NOT NULL)
            """)
        conn.commit()

    def add_go_term(self, upac, goid, name='', qualifier='', aspect=''):
        conn = self._get_conn()
        curs = conn.cursor()
        curs.execute("""
            INSERT INTO GoTerm(UniProtAccNum, GoId, Name, Qualifier, Aspect)
            VALUES (?, ?, ?, ?, ?)
            """, (upac, goid, name, qualifier, aspect))
        conn.commit()

    def delete_go_terms(self):
        conn = self._get_conn()
        curs = conn.cursor()
        curs.execute("""
            DELETE FROM GoTerm
            """)
        conn.commit()


    # Routines for PDB table

    def create_pdb_table(self):
        conn = self._get_conn()
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
        conn = self._get_conn()
        curs = conn.cursor()
        curs.execute("""
            VACUUM
            """)
        conn.commit()

    def dump_db(self):
        conn = self._get_conn()
        with codecs.open('icdb_dump.sql', 'w', 'utf-8') as dump_file:
            for line in conn.iterdump():
                dump_file.write('%s\n' % line)

    def print_db_stats(self):
        conn = self._get_conn()
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
        upac_list = self.get_protein_uniprot_accnums()
        for upac in upac_list:
            interaction_id_list = \
                self.get_interaction_ids_by_uniprot_accnum(upac)
            expression_id_list = \
                self.get_expression_ids_by_uniprot_accnum(upac)
            if interaction_id_list and expression_id_list:
                useful_data_count += 1
        print 'Number of proteins with interaction and expression = %d' \
            % useful_data_count
        upac_list = self.get_protein_uniprot_accnums()
        in_betse_count = 0
        in_betse_with_expr_count = 0
        in_betse_with_inter_count = 0
        for upac in upac_list:
            if self.is_in_betse(upac):
                in_betse_count += 1
                if self.exists_expression(upac):
                    in_betse_with_expr_count += 1
                if self.exists_interaction_by_uniprot_accnum(upac):
                    in_betse_with_inter_count += 1
        print 'Number of proteins in betse = %d' % in_betse_count
        print 'Number of proteins in betse with expression = %d' \
            % in_betse_with_expr_count
        print 'Number of proteins in betse with interaction = %d' \
            % in_betse_with_inter_count
