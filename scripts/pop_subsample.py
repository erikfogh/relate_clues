import argparse
import os
path_to_ids = "/faststorage/project/simons/data/1000Genomes/metainfo/"

parser = argparse.ArgumentParser(description = "script to make a single file containing all ids to remove")
parser.add_argument("ids_to_keep", nargs = "+", help = "files containing the ids to keep")
parser.add_argument("-o", help = "name and path to write the finished file")

args = parser.parse_args()

ids_to_remove = []
for f in os.listdir(path_to_ids):
    if f not in args.ids_to_keep and f.endswith("le.txt") and not f.startswith("all") and not f.startswith("."):
        ids_to_remove.append(f)

with open(args.o, "w") as outfile:
    for f in ids_to_remove:
        with open(os.path.join(path_to_ids, f)) as infile:
            outfile.write(infile.read())

print("Done with ", args.o)