import argparse
import os
import random
path_to_ids = "/faststorage/project/simons/data/1000Genomes/metainfo/"

parser = argparse.ArgumentParser(description="script to make a single file containing all ids to remove")
parser.add_argument("number", help="number of individuals to keep. Do not specify a number higher than the number of individuals")
parser.add_argument("-o", help="name and path to write the finished file")

args = parser.parse_args()

#I have just altered the code for pop_subsample, so this starting part isnt as direct as one could wish.
ids_to_remove = []
for f in os.listdir(path_to_ids):
    if f.endswith("le.txt") and not f.startswith("all") and not f.startswith("."):
        ids_to_remove.append(f)

with open(args.o+"_temp.txt", "w") as outfile:
    for f in ids_to_remove:
        with open(os.path.join(path_to_ids, f)) as infile:
            outfile.write(infile.read())


samples_left = int(args.number)
individuals_left = sum(1 for line in open(args.o+"_temp.txt"))
count = 0
print("We have {} individuals, and want to keep {}".format(individuals_left, samples_left))
f = open(args.o+".txt", "w").close()
with open(args.o+"_temp.txt", "r") as infile:
    for line in infile:
        prob = samples_left / individuals_left
        if random.random() <= prob:
            samples_left -= 1
            count += 1
        else:
            f = open(args.o+".txt", "a")
            f.write(line)
        individuals_left -= 1
os.remove(args.o+"_temp.txt")

print("Done with ", args.o)
print("Kept {} individuals".format(count))
