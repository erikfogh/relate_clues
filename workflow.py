from gwf import Workflow, AnonymousTarget
import os
import glob
from groups import Group

gwf = Workflow()

path_to_relate = "/home/eriks/baboondiversity/people/eriks/relate_clues_baboons/relate/relate_v1.1.7_x86_64_dynamic/"
genetic_map = "/home/eriks/baboondiversity/data/PG_panu3_recombination_map/genetic_map_chr{}.txt"
path_to_prepared = "steps/"
pop_information = {"all_individuals_prepared/": ["0.57e-8", "16000"]}
base_path = os.getcwd()
# "CEU_kept_prepared/": ["1.25e-8", "6000"],, "all_individuals_prepared/": ["1.25e-8", "30000"]  # Due to the larger number of samples, jobs fail with out-of-time or memory.
# for pop_information: key is name of folder in steps, values are mutation rate and effective population size
# chromosomes = list(range(1, 23))+["X"]
chromosomes = ["8"] # list(range(1, 23))+["X"]
job_name = "chr8"
l_d = []
for chrom in chromosomes:
    d = {}
    d["number"] = chrom
    d["genetic_map"] = genetic_map.format(chrom)
    l_d.append(d)

def full_relate(number, genetic_map, out_dir, pop_inf, prep_dir, base_path, relate_path):
    haps = base_path+"/"+prep_dir+"chrom{}.haps.gz".format(number)
    sample = base_path+"/"+prep_dir+"chrom{}.sample.gz".format(number)
    annot = base_path+"/"+prep_dir+"chrom{}.annot".format(number)
    dist = base_path+"/"+prep_dir+"chrom{}.dist.gz".format(number)
    m = pop_inf[0]
    n = pop_inf[1]
    o = "chrom{}".format(number)
    inputs = {"haps": haps, "sample": sample}
    Relate = os.path.join(relate_path, "bin/Relate")
    outputs = [out_dir+o+".anc"]
    options = {
        "cores": 7,
        "memory": "28g",
        "walltime": "24:00:00",
        "account": "baboondiversity"
    }
    spec = """
    cd {}
    {} --mode All -m {} -N {} --haps {} --sample {} --map {} -o {} --annot {} --dist {} --memory 20
    """.format(out_dir,
        Relate, m, n, haps, sample, genetic_map, o, annot,  dist)
    return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)


def estimate_pop_size(number, out_dir, pop_inf, prep_dir, relate_path):
    i = out_dir+"chrom{}.anc".format(number)
    o = out_dir+"chrom{}_popsize".format(number)
    m = pop_inf[0]
    inputs = [out_dir+"chrom{}.anc".format(number)]
    outputs = [o+".coal"]
    Relate = os.path.join(relate_path, "scripts/EstimatePopulationSize/EstimatePopulationSize.sh")
    poplabels = prep_dir+"chrom{}.poplabels".format(number)
    options = {
        "cores": 15,
        "memory": "60g",
        "walltime": "24:00:00",
        "account": "baboondiversity"
    }
    spec = """
    {} -i {} -m {} --poplabels {} -o {} --threshold 0 --num_iter 5 --years_per_gen 11 --threads 14
    """.format(Relate, i[:-4], m, poplabels, o)
    return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)


def detect_selection(number, in_dir, out_results, pop_inf, prep_dir, relate_path):
    i = in_dir+"chrom{}_popsize.coal".format(number)
    o = out_results+"chrom{}_selection".format(number)
    m = pop_inf[0]
    inputs = [i]
    outputs = [o+".freq"]
    Relate = os.path.join(relate_path, "scripts/DetectSelection/DetectSelection.sh")
    options = {
        "cores": 5,
        "memory": "20g",
        "walltime": "1:00:00",
        "account": "baboondiversity"
    }
    spec = """
    {} -i {} -m {} -o {}
    """.format(Relate, i[:-5], m, o)
    return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)


# def sample_branch_lengths(number, out_dir, pop_inf, prep_dir, relate_path):
#     i = out_dir+"chrom{}_popsize.coal".format(number)
#     o = out_dir+"chrom{}_selection".format(number)
#     m = pop_inf[0]
#     inputs = [i]
#     outputs = [o+".freq"]
#     Relate = os.path.join(relate_path, "scripts/DetectSelection/DetectSelection.sh")
#     options = {
#         "cores": 4,
#         "memory": "2",
#         "walltime": "1:00:00"
#     }
#     spec = """
#     {} -i {} -m {} -o {}
#     """.format(Relate, i[:-4], m, o)
#     print(spec)
#     return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)


for pop in pop_information:
    pop_name = pop[:-10]
    out_dir = "steps/{}_relate/".format(pop_name)
    out_results = "results/{}_relate/".format(pop_name)
    os.makedirs(out_dir, exist_ok=True)
    os.makedirs(out_results, exist_ok=True)
    with Group(gwf, suffix=pop_name) as g:
        relate = g.map(full_relate, l_d, name=job_name+"relate", extra={
            "out_dir": out_dir, "pop_inf": pop_information[pop],
            "prep_dir": path_to_prepared+pop, "base_path": base_path, "relate_path": path_to_relate
        })
        popsize = g.map(estimate_pop_size, chromosomes, name=job_name+"popsize", extra={
            "out_dir": out_dir, "pop_inf": pop_information[pop],
            "prep_dir": path_to_prepared+pop, "relate_path": path_to_relate
        })
        g.map(detect_selection, chromosomes, name=job_name+"detect_selection", extra={
            "in_dir": out_dir, "out_results": out_results, "pop_inf": pop_information[pop],
            "prep_dir": path_to_prepared+pop, "relate_path": path_to_relate
        })

#
# test_input = gwf.target_from_template(
#     name="CEU_chrom22_test",
#     template=full_relate(
#         m="1.25e-8", n="30000", prep_dir="steps/CEU_kept_prepared/", number="22",
#         genetic_map="data/recombination_maps/genetic_map_chr22_combined_b37.txt",
#         relate_path=path_to_relate, o="steps/CEU_chrom22_test"
#     )
# )

# test_popsize = gwf.target_from_template(
#     name="CEU_chrom22_popsize",
#     template=estimate_pop_size(
#         m="1.25e-8", i=test_input.outputs[0],
#         prep_dir="steps/CEU_kept_prepared/", number="22",
#         relate_path=path_to_relate, o="CEU_chrom22_poptest"
#     )
# )

# test_detection = gwf.target_from_template(
#     name="CEU_chrom22_selection",
#     template=detect_selection(
#         m="1.25e-8", i=test_input.outputs[0], c=test_popsize.outputs[0],
#         relate_path=path_to_relate, o="CEU_chrom22_seletest"
#     )
# )
