#! /usr/bin/python

import sys
import FileUtility
import Error
import sqlite3


class Genomes_to_taxon():
    def __init__(self):
        self._converter = sqlite3.connect('/rap/nne-790-ab/projects/pplante/Git_clone/laughing-nemesis/GenomesToTaxon.sqlite')
        self._cursor = self.converter.cursor()

    @property
    def converter(self):
        return self._converter

    @property
    def cursor(self):
        return self._cursor

    def populate(self, file_name1, taxon_name):
        # Create table...
        curs = self.cursor
        curs.execute("create table GENOME_TO_TAXON (GENOMES integer, TAXON integer, PRIMARY KEY (GENOMES))")
        self.converter.commit()
        # Check is file is ok
        FileUtility.isValid(file_name1)
        sys.stderr.write("counting lines to prepare converter\n")
        n_lines = FileUtility.countLines(file_name1)
        sys.stderr.write(str(n_lines)+" to read\n")
        readed = 0
        # populate the DB
        for line in open(file_name1):
            genome = int(line.split()[0])
            taxon = int(line.split()[1])
            curs.execute('insert into GENOME_TO_TAXON (GENOMES, TAXON) values (?, ?)', (genome, taxon))
            readed += 1
            if readed % 10000 == 0:
                self.converter.commit()
                sys.stderr.write(str(readed)+" lines readed out of "+str(n_lines)+"\n")
        FileUtility.isValid(taxon_name)
        sys.stderr.write("counting lines to prepare converter\n")
        n_lines = FileUtility.countLines(taxon_name)
        sys.stderr.write(str(n_lines)+" to read\n")
        readed = 0
        curs.execute("create table TAXON_NAMES (TAXON integer, NAME VARCHAR(100), RANK VARCHAR(30), "
                     "PRIMARY KEY (TAXON))")
        for line in open(taxon_name):
            taxon = int(line.split("\t")[0])
            name = str(line.split("\t")[1])
            rank = str(line.split("\t")[2]).rstrip("\n")
            curs.execute('insert into TAXON_NAMES (TAXON, NAME, RANK) values (?, ?, ?)', (taxon, name, rank))
            readed += 1
            if readed % 10000 == 0:
                self.converter.commit()
                sys.stderr.write(str(readed)+" lines readed out of "+str(n_lines)+"\n")

    def get_taxon_name(self, taxon_id):
        curs = self.cursor
        curs.execute("select * from TAXON_NAMES where TAXON = ?", (taxon_id,))
        rows = curs.fetchall()
        if len(rows) != 1:
            raise ValueError("Taxon "+str(taxon_id)+" is not in the converter")
        else:
            return rows[0][1]

    def get_taxon_rank(self, taxon_id):
        curs = self.cursor
        curs.execute("select * from TAXON_NAMES where TAXON = ?", (taxon_id,))
        rows = curs.fetchall()
        if len(rows) != 1:
            raise ValueError("Taxon "+str(taxon_id)+" is not in the converter")
        else:
            return rows[0][2]



    def close(self):
        self.converter.close()

    def convert_to_taxon(self, genome):
        curs = self.cursor
        curs.execute("select * from GENOME_TO_TAXON where GENOMES = ?", (genome,))
        rows = curs.fetchall()
        if len(rows) != 1:
            raise ValueError("Genome "+str(genome)+" is not in the converter")
        else:
            return rows[0][1]

    def genome_is_valid(self, genome):
        curs = self.cursor
        curs.execute("select * from GENOME_TO_TAXON where GENOMES = ?", (genome,))
        rows = curs.fetchall()
        if len(rows) <= 0:
            return False
        return True


#test main      
if __name__=="__main__":
    if len(sys.argv) != 3:
        sys.stderr.write("To prepare the database, you must provide the Genomes_to_taxon.tsv file and the "
                         "Taxon-names.tsv file\n")
        sys.stderr.write("Ex: GenomesToTaxon.py Genomes_to_taxon.tsv Taxon-names.tsv file\n")
        sys.exit(0)
    conv = Genomes_to_taxon()
    conv.populate(sys.argv[1], sys.argv[2])
    conv.converter.close()

