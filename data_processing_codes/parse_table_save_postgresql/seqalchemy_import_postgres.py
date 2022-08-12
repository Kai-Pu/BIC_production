import pandas as pd
import sys
from sqlalchemy import create_engine

## import DB setting 
DATABASE = sys.argv[1]
USER = sys.argv[2]
PASSWORD = sys.argv[3]
ADDRESS = sys.argv[4]
PORT = sys.argv[5]
DB_NAME = sys.argv[6]

## import data and DB table setting
DATA_PATH = sys.argv[7]
DB_TABLE_NAME = sys.argv[8]


# dialect[+driver]://user:password@host:port/dbname
#URL = 'postgresql://[USER]:[PASSWORD]@[Address]:[PORT]/[DBNAME]'
URL = "{}://{}:{}@{}:{}/{}".format(DATABASE, USER, PASSWORD, ADDRESS, PORT, DB_NAME)
print(URL)

engine = create_engine(URL, echo=True)
data = pd.read_csv(DATA_PATH, sep="\t", encoding='utf-8')
data.to_sql(DB_TABLE_NAME, con=engine, index=False, if_exists="append")
engine.dispose()


print('Done!')