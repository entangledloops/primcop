#!/usr/bin/python

sql = """SELECT
	CAST ('max_papa_score','max_prima_score', 'fold_index_score' AS DOUBLE PRECISION);
	"""

import psycopg2
from config import config

def insert_sequences():

	conn = None
	try:
		# read the connection parameters
		params = config()
		# connect to the PostgreSQL server
		conn = psycopg2.connect(**params)
		cur = conn.cursor()
		# open csv file to read data into table
		fin1 = open('/home/wrivera/primcop/sequences.csv', 'r')
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
	insert_sequences()