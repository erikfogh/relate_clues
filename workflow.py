from gwf import Workflow, AnonymousTarget
import os
import glob
from groups import Group

gwf = Workflow()

path_to_relate = "/faststorage/home/eriks/relate-clues/relate_v1.1.2_x86_64_dynamic/"
pop_files = "data/pops/"
use_all = True
pop_list = glob.glob(pop_files + "*.txt")
path_to_prepared = "steps/"

prep_list = {"CEU": "steps/CEU_kept_prepared"}


def full_relate(m, n, prep_dir, number, genetic_map, relate_path, o):
    haps = prep_dir+"chrom{}.haps.gz".format(number)
    sample = prep_dir+"chrom{}.sample.gz".format(number)
    annot = prep_dir+"chrom{}.annot".format(number)
    dist = prep_dir+"chrom{}.dist.gz".format(number)
    inputs = {"haps": haps, "sample": sample}
    Relate = os.path.join(relate_path, "bin/Relate")
    outputs = [o+".anc"]
    options = {
        "cores": 8,
        "memory": "20g",
        "walltime": "2:00:00"
    }
    spec = """
    {} --mode All -m {} -N {} --haps {} --sample {} --map {} -o {} --annot {} --dist {} --memory 17
    """.format(Relate, m, n, haps, sample, genetic_map, o, annot,  dist)
    return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)


def estimate_pop_size(m, i, prep_dir, number, o, relate_path):
    inputs = [i]
    outputs = [o+".coal"]
    Relate = os.path.join(relate_path, "scripts/EstimatePopulationSize/EstimatePopulationSize.sh")
    poplabels = prep_dir+"chrom{}.poplabels".format(number)
    options = {
        "cores": 10,
        "memory": "20g",
        "walltime": "10:00:00"
    }
    spec = """
    {} -i {} -m {} --poplabels {} -o {} --threshold 0 --num_iter 5
    """.format(Relate, i[:-4], m, poplabels, o)
    print(spec)
    return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)


def detect_selection(m, i, o, relate_path):
    inputs = [i]
    outputs = [o+".freq"]
    Relate = os.path.join(relate_path, "/scripts/DetectSelection/DetectSelection.sh")
    options = {
        "cores": 10,
        "memory": "20g",
        "walltime": "10:00:00"
    }
    spec = """
    {} -i {} -m {} -o {}
    """.format(Relate, i[:-4], m, o)
    print(spec)
    return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)


test_input = gwf.target_from_template(
    name="CEU_chrom22_test",
    template=full_relate(
        m="1.25e-8", n="30000", prep_dir="steps/CEU_kept_prepared/", number="22",
        genetic_map="data/recombination_maps/genetic_map_chr22_combined_b37.txt",
        relate_path=path_to_relate, o="CEU_chrom22_test"
    )
)

test_popsize = gwf.target_from_template(
    name="CEU_chrom22_popsize",
    template=estimate_pop_size(
        m="1.25e-8", i=test_input.outputs[0],
        prep_dir="steps/CEU_kept_prepared/", number="22",
        relate_path=path_to_relate, o="CEU_chrom22_poptest"
    )
)
