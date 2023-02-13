from gwf import Workflow, AnonymousTarget
import os
from groups import Group

gwf = Workflow()

# 

path_to_relate = "/home/eriks/baboondiversity/people/eriks/second_analysis_baboons/relate/relate_v1.1.7_x86_64_dynamic/"
in_dir_prep = "steps/all_individuals_prepared/"
in_dir_relate = "steps/all_individuals_relate/"
poplabels = "data/pops/all_females_8cluster.sample"
out_path = "results/all_females/"
chromosome = "X"
job_name = "chrX"


def plot_trees(loc, name, in_dir_prep, in_dir_relate, poplabels, out_path, relate_path):
    anc = in_dir_relate+name+"_popsize.anc.gz"
    mut = in_dir_relate+name+"_popsize.mut.gz"
    haps = in_dir_prep+name+".haps.gz"
    sample = in_dir_prep+name+".sample.gz"
    inputs = [anc, haps]
    outname = name+"_"+str(loc)
    outputs = [out_path+outname+".pdf"]
    Relate = os.path.join(relate_path, "scripts/TreeView/TreeView.sh")
    options = {
        "cores": 4,
        "memory": "30g",
        "walltime": "1:00:00",
        "account": "baboondiversity"
    }
    spec = """
    {} --anc {} --mut {} --haps {} --sample {} --poplabels {} \
     --bp_of_interest {}  --years_per_gen 11 -o {}
    """.format(Relate, anc, mut, haps, sample, poplabels,
               loc, out_path+outname)
    return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)


locs = [137576566,
 135090332,
 137266109,
 135282236,
 135805937,
 136164193,
 135659139,
 136717983,
 136895491,
 136800825,
 136387268,
 136822194,
 136092096,
 136082477,
 136280537,
 138160118,
 136098298,
 135090234,
 137558239,
 136441507]# [17500000, 18000000, 18400022, 19000000] #[16623680, 16607974, 17000000, 66182731, 74015344, 110716419]


os.makedirs(out_path, exist_ok=True)
with Group(gwf, suffix=job_name) as g:
    g.map(plot_trees, locs, name="plot_trees", extra={
            "name": "chrom"+chromosome, "in_dir_prep": in_dir_prep, "in_dir_relate": in_dir_relate,
            "poplabels": poplabels, "out_path": out_path, "relate_path": path_to_relate
        })
