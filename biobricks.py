import biobricks as bb
import pandas as pd
import pyspark


brick = bb.assets('bindingdb')
print(brick)
brick_dict = vars(brick)

spark = pyspark.sql.SparkSession.builder.appName("cpdat_analysis").config("spark.sql.debug.maxToStringFields", 1000).getOrCreate()
# spark.sql.debug.maxToStringFields()
for key in brick_dict:
    print(key, "->", brick_dict[key], "\n")
    df = spark.read.parquet(brick_dict[key])
    df.show()

print("------------------------------------------------------------")
print("------------------------------------------------------------")
#
# import sqlite3
#
# chembl_brick = bb.assets('chembl')
# chembl = vars(chembl_brick)
#
# for key in chembl:
#     print(key, "->", chembl[key], "\n")
#     sqlite_path = chembl[key]
#     conn = sqlite3.connect(sqlite_path)
#     cursor = conn.cursor()
#
#     # List all tables
#     cursor.execute("SELECT name FROM sqlite_master WHERE type='table';")
#     tables = cursor.fetchall()
#     print("Tables:", tables)
