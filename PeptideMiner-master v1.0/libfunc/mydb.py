import re
import sqlite3

"""
Interface to SQLite DB. 
"""

class mysql:
	def __init__ (self):

		self.con = sqlite3.connect('sqlite.db')
		self.cur = self.con.cursor()
		self.test = False
		self.verbose = False
		self.e = None

	def last_error (self):
		print(f"Last error: {self.e}")

	def execute (self,query):
		try:
			self.cur.execute(query)
			self.con.commit()
			return 1
		except sqlite3.Error as e:
			self.e = e
			print(f"ERROR MYSQL> {e}")
			print(f" << the query was : {query} >>")
			return 0

	def __build_query__ (self,fields='',tables='',where='',other=''):
		if type(fields) != str:
			fields = ",".join(fields)

		query = f"select {fields}")
		if type(tables) != str:
			if type(tables) == dict:
				listtab = []
				for (k,v) in tables.iteritems():
					listtab.append(f"{k} {v}")
				tables = ",".join(listtab)
			else:
				tables = ",".join(tables)

		if tables != '':
			query = f"{query} from {tables}"

		if type(where) != str:
			if type(where) == dict:
				listw = []
				for (k,v) in where.iteritems():
					if re.search('^[a-z]{1,4}\.[A-za-z_]+$',str(v)): # refer another field
						listw.append(f"{k}={v}")
					else:
						listw.append("{k}=\"{1}\"".format(k,re.sub("\"","\\\"",str(v))))
				where = " and ".join(listw)
			else:
				where = " and ".join(where)
		if where != '':
			query = f"{query} where {where}"
		if other != '':
			query = f"{query} {other}"
		if self.verbose:
			print(query)
		return query

	def onerow (self,fields='',tables='',where='',other=''):
		self.execute(self.__build_query__(fields,tables,where,other))
		return self.cur.fetchone()
	
	def nextrow (self):
		return self.cur.fetchone()

	def onevalue (self,fields='',tables='',where='',other=''):
		row = self.onerow(fields,tables,where,other)
		if row:
			return row[0]
		return None

	def __row_dict__ (self,fields,res):
		if not res:
			return None
		resd = dict()
		if type(fields) == str:
			fields = re.split(" *, *",fields)
		for i in range(0,len(fields)):
			resd[fields[i]] = res[i]
		return resd

	def onerow_dict (self,fields='',tables='',where='',other=''):
		res = self.onerow(fields,tables,where,other)
		if not res:
			return None
		return self.__row_dict__(fields,res)

	def nextrow_dict (self,fields):
		return self.__row_dict__(fields,self.cur.fetchone())
	
	def selectall (self,fields='',tables='',where='',other=''):
		if self.execute(self.__build_query__(fields,tables,where,other)):
			return self.cur.fetchall()
		return None

	def selectall_dict (self,fields='',tables='',where='',other=''):
		res = self.selectall(fields,tables,where,other)
		if not res:
			return None
		resd = []
		if type(fields) == str:
			fields = re.split(" *, *",fields)
		for r in res:
			rd = dict()
			for i in range(0,len(fields)):
				rd[fields[i]] = r[i]
			resd.append(rd)
		return resd
	
	def insert (self,table,values):
		entered_values = dict()
		for (k,v) in values.iteritems():
			if v == None:
				continue
			elif not re.search("^(LAST_INSERT_ID|NOW)\(\)",str(v)):
				v = re.sub('\\\\','\\\\\\\\',str(v))
				values[k] = "'{0}'".format(re.sub('\'','\\\\\'',str(v)))
			entered_values[k] = values[k]

		sqlIns  = f"insert into {table} "
		sqlIns += f" ({','.join(entered_values.keys())}) "
		sqlIns += f" values ({','.join([str(v) for v in entered_values.values()])})"
		if not self.test:
			self.execute(sqlIns)

		if self.test and self.verbose:
			print(sqlIns)
			

	def delete (self,table,where=''):
		query = "delete from {0}".format(table)
		if where != '':
			query = "{0} where {1}".format(query,where)
		if not self.test:
			self.execute(query)
		if self.test and self.verbose:
			print(query)
	
	def update (self,table,values,where=''):
		query = "update {0} set ".format(table)
		updates = []
		for c in values.keys():
			if type(values[c]) != str:
				m = False
			else:
				m = re.search("ff(.*)ff",values[c]) #code for a field
			if m:
				updates.append("{0}={1}".format(c,m.group(1)))
			else:
				val = re.sub("([^\\\])'","\\1\\'",str(values[c]))
				updates.append("{0}='{1}'".format(c,val))
		query = "update {0} set {1}".format(table,",".join(updates))
		if where:
			query = "{0} where {1}".format(query,where)
		if self.verbose:
			print(query)
		if not self.test:
			self.execute(query)
	
	def alter (self,table,command):
		query = "alter table {0} {1}".format(table,command)
		if not self.test:
			self.execute(query)
		if self.test and self.verbose:
			print(query)

	def findminid (self,table,field):
		maxi = 1
		while self.onevalue(field,table,"{0}={1}".format(field,maxi)) != None:
			maxi = maxi+1
		return maxi

	def setminautoincrement (self,table,field):
		self.alter(table,f"AUTO_INCREMENT={self.findminid(table,field)}")

	def create_if_not (self,db):
		query = f"create database if not exists {db}"
		if not self.test:
			self.execute(query)
	
		if self.test and self.verbose:
			print(query)

	def use (self,db):
		query = f"use {db}"
		if not self.test:
			self.execute(query)

		if self.test and self.verbose:
			print(query)
