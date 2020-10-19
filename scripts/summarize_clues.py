import argparse
import pandas as pd

parser = argparse.ArgumentParser(description="gather all the results created by clues")
parser.add_argument("-i", nargs="+", help="list of files to summarize")
parser.add_argument("-o", help="name and path to write the finished file")

args = parser.parse_args()

count = 0
for i in args.i:
    count += 1

print(count)
