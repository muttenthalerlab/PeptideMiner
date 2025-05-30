import os,re
import sqlite3

import logging
logger = logging.getLogger(__name__)

"""
Interface to SQLite DB. 
"""

#=======================================================================
class sql_connector():
#=======================================================================	
    def __init__ (self, data_file= 'sqlite.db', sql_create= 'sqlite_create.sql'):
        self.database = data_file

        new_db = not os.path.isfile(self.database)

        self.con = sqlite3.connect(self.database)
        self.cur = self.con.cursor()

        # If new database file - create 
        if new_db:
            if os.path.isfile(sql_create):
                logger.info(f" [SQLite] New Database {self.database} using {sql_create}")
                with open(sql_create) as f:
                    sql_script = f.read()
                self.cur.executescript(sql_script)
        else:
            logger.info(f" [SQLite] Database {self.database} ")

        self.test = False
        self.verbose = False
        self.e = None

    def last_error (self):
        logger.error(f" [SQLite] Last Error: {self.e}")

    def execute (self,query):
        try:
            self.cur.execute(query)
            self.con.commit()
            return 1
        except sqlite3.Error as e:
            self.e = e
            logger.error(f" [SQLite] Error: {e} from query {query}")
            return 0

    def __build_query__ (self,tables='',fields='',where='',other=''):
        if type(fields) != str:
            fields = ",".join(fields)

        query = f"select {fields}"

        if type(tables) != str:
            if type(tables) == dict:
                listtab = []
                for k in tables:
                    listtab.append(f"{k} {tables[k]}")
                tables = ",".join(listtab)
            else:
                tables = ",".join(tables)
        if tables != '':
            query = f"{query} from {tables}"

        if type(where) != str:
            if type(where) == dict:
                listw = []
                for k in where:
                    if re.search('^[a-z]{1,4}\.[A-za-z_]+$',str(where[k])): # refer another field
                        listw.append(f"{k}={where[k]}")
                    else:
                        listw.append("{0}=\'{1}\'".format(k,re.sub("\"","\\\"",str(where[k]))))
                where = " and ".join(listw)
            else:
                where = " and ".join(where)
                
        if where != '':
            query = f"{query} where {where}"
        if other != '':
            query = f"{query} {other}"
        if self.verbose:
            logger.info(query)
        return query

    def onerow (self,tables='',fields='',where='',other=''):
        self.execute(self.__build_query__(tables,fields,where,other))
        return self.cur.fetchone()
	
    def nextrow (self):
        return self.cur.fetchone()

    def onevalue (self,tables='',fields='',where='',other=''):
        row = self.onerow(tables,fields,where,other)
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

    def onerow_dict (self,tables='',fields='',where='',other=''):
        res = self.onerow(tables,fields,where,other)
        if not res:
            return None
        return self.__row_dict__(fields,res)

    def nextrow_dict (self,fields):
        return self.__row_dict__(fields,self.cur.fetchone())
	
    def selectall (self,tables='',fields='',where='',other=''):
        if self.execute(self.__build_query__(tables,fields,where,other)):
            return self.cur.fetchall()
        return None

    def selectall_dict (self,tables='',fields='',where='',other=''):
        res = self.selectall(tables,fields,where,other)
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
        for k in values:
            if values[k] == None:
                continue
            elif not re.search("^(LAST_INSERT_ID|NOW)\(\)",str(values[k])):
                values[k] = re.sub('\\\\','\\\\\\\\',str(values[k]))
                values[k] = "'{0}'".format(re.sub('\'','\\\\\'',str(values[k])))
            entered_values[k] = values[k]

        sqlIns  = f"insert into {table} "
        sqlIns += f" ({','.join(entered_values.keys())}) "
        sqlIns += f" values ({','.join([str(v) for v in entered_values.values()])})"
        #print(f" [-----] {sqlIns}")
        if not self.test:
            self.execute(sqlIns)

        if self.test and self.verbose:
            logger.info(sqlIns)
			

    def delete (self,table,where=''):
        query = "delete from {0}".format(table)
        if where != '':
            query = "{0} where {1}".format(query,where)
        if not self.test:
            self.execute(query)
        if self.test and self.verbose:
            logger.info(query)
	
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
            logger.info(query)
        if not self.test:
            self.execute(query)
	
    def alter (self,table,command):
        query = "alter table {0} {1}".format(table,command)
        if not self.test:
            self.execute(query)
        if self.test and self.verbose:
            logger.info(query)

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
            logger.info(query)

    def use (self,db):
        query = f"use {db}"
        if not self.test:
            self.execute(query)

        if self.test and self.verbose:
            logger.info(query)