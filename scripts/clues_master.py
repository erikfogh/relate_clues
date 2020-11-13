import os
import argparse
import pandas as pd
from multiprocessing import Pool
import re
import glob
import subprocess

# I have chosen to only have the iterable defined as a local variable in the function call,
# while all other arguments are taken either from the temp table or global variables.

parser = argparse.ArgumentParser(description="script to handle multiple clues jobs, and output the results")
parser.add_argument("chunk", type=int, help="which chunk this job should handle")
parser.add_argument("total_chunks", type=int, help="number of chunks")
parser.add_argument("-o", help="path to write the finished file")
parser.add_argument("-i", help="name and path to the relate files without suffix")
parser.add_argument("-m", help="mutation rate")
parser.add_argument("-r", help="path to relate")
parser.add_argument("-b", help="time bin")
parser.add_argument("--burnin", help="burnin for clues", default=50)
parser.add_argument("--thin", help="thinning for clues", default=1)
parser.add_argument("--mcmc", help="hos many trees relate should run the mcmc for", default=150)

args = parser.parse_args()
temp_table = pd.read_table(args.o+"clues_table_temp.txt", sep='\s+')
table_length = len(temp_table)
start = -(-(args.chunk-1)*table_length//args.total_chunks)
end = -(-args.chunk*table_length//args.total_chunks)
number_of_cores = int(os.environ['SLURM_CPUS_PER_TASK'])
snps = temp_table.iloc[start:end]["pos"].to_list()
chunk_output = args.o+"tmp/chunk{}_*.txt".format(args.chunk)

print("There are {} cores".format(number_of_cores))


def relate_clues(snp):
    tempchunk = args.o+"tmp/chunk{}_snp{}".format(args.chunk, snp)
    coal = args.i+".coal"
    dist = args.i+".dist"
    daf = float(temp_table.loc[temp_table["pos"] == snp, "daf"])
    spec1 = """
    {}bin/RelateExtract --mode AncMutForSubregion  --first_bp {} --last_bp {} --anc {}.anc --mut {}.mut -o {}
    """.format(args.r, snp, snp, args.i, args.i, tempchunk)
    subprocess.run(args=spec1, shell=True)
    spec2 = """
    {}/bin/RelateCoalescentRate --mode SampleBranchLengths -m {} --coal {} --num_samples {} -i {} -o {} --format {} --dist {}
    """.format(args.r, args.m, coal, args.mcmc, tempchunk, tempchunk, "b", dist)
    subprocess.run(args=spec2, shell=True)
    spec3 = """
    python inference.py --times {} --coal {} --timeBins {} --popFreq {} --burnin {} --out {} > {}_output.txt
    """.format(tempchunk, coal, args.b, daf, args.burnin, tempchunk, tempchunk)
    subprocess.run(args=spec3, shell=True)
    print("Finished snp {}".format(snp))


with Pool(number_of_cores, maxtasksperchild=5) as pool:
    results = pool.map(relate_clues, snps)
    pool.close()
    pool.join()
tempfiles = glob.glob(chunk_output)
snps = []
clues_LR = []
for path in tempfiles:
    snp = re.match(".+snp(\d+)",path)
    snps.append(int(snp.group(1)))
    f = open(path, "r").read()
    logLR = re.search("logLR: (\d+.\d+)", f)
    clues_LR.append(float(logLR.group(1)))
clues_data = pd.DataFrame(list(zip(snps, clues_LR)), columns = ["pos", "clues_LR"])
clues_data.to_csv(args.o+"chunk{}_table.txt".format(args.chunk), index=False, sep=" ")
print("All jobs are done")

# {}bin/RelateExtract --mode AncMutForSubregion  --first_bp {} --last_bp {} --anc {}.anc --mut {}.mut -o {}
#     {}/bin/RelateCoalescentRate --mode SampleBranchLengths -m {} --coal {} --num_samples {} -i {} -o {} --format {} --dist {}
#     cd clues/
#     python inference.py --times {} --coal {} --out {} --timeBins {} > {}_output.txt
# .format(relate_path, snp, snp, i[:-5], i[:-5], tmp,
#                relate_path, m, i, 150, tmp, tmp, "b", dist,
#                tmp, coal, o, timebins, o)
