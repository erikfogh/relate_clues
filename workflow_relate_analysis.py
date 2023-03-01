from gwf import Workflow, AnonymousTarget
import os
import pandas as pd
from groups import Group

gwf = Workflow()
poplabels_path = "data/pops/all_females_8cluster.sample"
poplabels = pd.read_csv(poplabels_path, sep=" ",
                        names=["i", "ID", "POP", "GROUP", "SEX"], header=0)
tree_path = "steps/all_individuals_relate/chromX_tskit.trees"
out_path = "steps/relate_coal_ordering/"
chromosomes = ["X"]  # list(range(1, 23))+["X"]
job_name = "chrX"
ID_list = poplabels.ID
l_d = []


def relate_coal_ordering(ID, tree_path, poplabels_path, out_path):
    inputs = [tree_path]
    outputs = [out_path + ID + ".txt"]
    options = {
        "cores": 2,
        "memory": "14g",
        "walltime": "4:00:00",
        "account": "baboondiversity"
    }
    spec = """
    python scripts/relate_coal_ordering.py -t {} -i {} -p {} -o {}
    """.format(tree_path, ID, poplabels_path, out_path)
    print(spec)
    return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)


os.makedirs(out_path, exist_ok=True)

gwf.map(relate_coal_ordering, ID_list, extra={"tree_path": tree_path, "poplabels_path": poplabels_path,
                                              "out_path": out_path})