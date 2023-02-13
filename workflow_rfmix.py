from gwf import Workflow, AnonymousTarget
import os
from groups import Group

gwf = Workflow()

### Input paths
# Uses the baboondiversity_v2 environment. Requirements are mostly bcftools and rfmix.

path_to_input = "/home/eriks/baboondiversity/data/PG_panu3_phased_chromosomes_4_7_2021/rfmix_input/"
query_file = path_to_input + "all_chr_females_query.bcf"
ref_file = path_to_input + "all_chr_females_reference.bcf"
sample_map = path_to_input + "sample_map_females.txt"
genetic_map = path_to_input + "all_chr_recombination_map.txt"
path_to_output = "steps/rfmix_females_gen50/"
base_path = os.getcwd()
chromosomes = list(range(1, 21)) + ["X"]
# list(range(1, 21))+

### Gwf functions

def rfmix(chrom, query, reference, sample_map, genetic_map, output_path):
    output = output_path + "chr" + str(chrom)
    inputs = [query, reference, sample_map, genetic_map]
    outputs = [output + ".msp.tsv"]
    options = {
        "cores": 10,
        "memory": "200g",
        "walltime": "16:00:00",
        "account": "baboondiversity"
    }
    spec = """
    rfmix -f {} -r {} -m {} -g {} -o {} --chromosome=chr{} -e 3 -G 50 --reanalyze-reference
    """.format(query, reference, sample_map, genetic_map, output, chrom)
    return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)


### Function calls

os.makedirs(path_to_output, exist_ok=True)

gwf.map(rfmix, chromosomes, name="rfmix_all",
                extra={"query": query_file, "reference": ref_file,
                       "sample_map": sample_map, "genetic_map": genetic_map,
                       "output_path": path_to_output})


