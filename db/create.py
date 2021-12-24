#!/usr/bin/python

sql = """SELECT
	CAST ('fold_index_score',max_papa_score','max_prima_score' AS DOUBLE PRECISION);
	"""

import psycopg2
from config import config

def create_and_pop_tables():
	""" create  and populate tables with .csv data in the PostgreSQL database"""
	commands = (
		"""
		CREATE TABLE contributors (
				institution_name VARCHAR(100),
				institution_id SERIAL NOT NULL,
				affiliated_lab VARCHAR(100),
				lab_id SERIAL NOT NULL PRIMARY KEY
		)
		""",
		"""
		CREATE TABLE sequences (
				sequence_id VARCHAR(25),
				sequence_number SERIAL NOT NULL PRIMARY KEY,
				seq_repeat VARCHAR(25),
				max_papa_score VARCHAR (15),
				max_prima_score VARCHAR(15),
				aa_sequence VARCHAR(1000)
		)
		""",
		"""
		CREATE TABLE contributor_sequences (
				lab_id INTEGER NOT NULL,
				sequence_number INTEGER NOT NULL,
				PRIMARY KEY (lab_id,sequence_number),
				FOREIGN KEY (lab_id)
					REFERENCES contributors (lab_id)
					ON UPDATE CASCADE ON DELETE CASCADE,
				FOREIGN KEY (sequence_number)
					REFERENCES contributors (lab_id)
					ON UPDATE CASCADE ON DELETE CASCADE
		)
		""")
	conn = None
	try:
		# read the connection parameters
		params = config()
		# connect to the PostgreSQL server
		conn = psycopg2.connect(**params)
		cur = conn.cursor()
		# create table one by one
		for command in commands:
			cur.execute(command)
		
		# open csv file to read data into table
		fin0 = open('/home/ubuntu/Documents/primcop.csv', 'r')
		#skip header
		next(fin0)
		# copy table data in from csv file
		cur.copy_from(fin0, 'contributors', columns=('institution_name','affiliated_lab'), sep=",")
		#close csv file
		fin0.close()
		# open csv file to read data into table
		fin1 = open('/home/ubuntu/Documents/hwseq.csv', 'r')
		#skip header
		next(fin1)
		# copy table data in from csv file
		cur.copy_from(fin1, 'sequences', null ='\\N', columns=('sequence_id', 'seq_repeat', 'max_papa_score', 'max_prima_score','aa_sequence'), sep=",")
		#close csv file
		fin1.close()
		# commit the changes
		conn.commit()
		# close communication with the PostgreSQL database server
		cur.close()

	except (Exception, psycopg2.DatabaseError) as error:
		print(error)
	finally:
		if conn is not None:
			conn.close()


if __name__ == '__main__':
	create_and_pop_tables()