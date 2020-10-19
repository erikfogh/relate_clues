import os
import argparse


parser = argparse.ArgumentParser(description="script to handle multiple clues jobs, and output the results")
parser.add_argument("chunk", type=int, help="which chunk this job should handle")
parser.add_argument("total_chunks", type=int, help="number of chunks")
parser.add_argument("total_jobs", type=int, help="number of test to solve")
parser.add_argument("-o", help="name and path to write the finished file")


# {}bin/RelateExtract --mode AncMutForSubregion  --first_bp {} --last_bp {} --anc {}.anc --mut {}.mut -o {}
#     {}/bin/RelateCoalescentRate --mode SampleBranchLengths -m {} --coal {} --num_samples {} -i {} -o {} --format {} --dist {}
#     cd clues/
#     python inference.py --times {} --coal {} --out {} --timeBins {} > {}_output.txt
# .format(relate_path, snp, snp, i[:-5], i[:-5], tmp,
#                relate_path, m, i, 150, tmp, tmp, "b", dist,
#                tmp, coal, o, timebins, o)
