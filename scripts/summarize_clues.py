import argparse
import pandas as pd
import glob

parser = argparse.ArgumentParser(description="gather all the results created by clues")
parser.add_argument("-i", help="input dir")
parser.add_argument("-o", help="name and path to write the finished file")

args = parser.parse_args()

tempfiles = glob.glob(args.i+"chunk*")
df_list = []
for path in tempfiles:
    f = pd.read_table(path, sep=" ")
    df_list.append(f)

df_clues = pd.concat(df_list)
df_relate = pd.read_table(args.i+"clues_table_temp.txt", sep=" ")

df_full = pd.merge(df_relate, df_clues)
print("Is there any NaN values?: ", df_full.isnull().values.any())
df_full.to_csv(args.o, index=False, sep=" ")
