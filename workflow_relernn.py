from gwf import Workflow, AnonymousTarget
import os
import glob
from groups import Group

gwf = Workflow()

### Input paths
# Uses the recomb environment.

path_to_input = "/home/eriks/baboondiversity/data/PG_panu3_phased_chromosomes_4_7_2021/relernn_input/"
path_to_output = "steps/relernn/"
genome_sizes = path_to_input + "autosome_sizes.bed"
mask_bed = path_to_input + "autosomes.sorted.bed"
mu = "0.57e-8"
base_path = os.getcwd()
#chromosomes = list(range(1, 21))

population_l = ["Dendro", "Arusha", "Ngorongoro"]

### Gwf functions

def relernn(name, input_dir, output_dir, GENOME, MASK, MU):
    project_dir = output_dir + name + "/"
    output = output_dir + name + "/" + name + ".PREDICT.txt"
    VCF = input_dir + name + ".vcf"
    inputs = [VCF, GENOME, MASK]
    outputs = [output]
    options = {
        "cores": 4,
        "memory": "30g",
        "walltime": "10:00:00",
        "account": "baboondiversity"
    }
    spec = """
    ReLERNN_SIMULATE --vcf {VCF} --genome {GENOME} \
    --mask {MASK} --projectDir {project_dir} --assumedMu {MU} \
    --nTrain 20000 --nVali 1000 --nTest 1000
    ReLERNN_TRAIN --projectDir {project_dir} --nEpochs 2 --nValSteps 2
    ReLERNN_PREDICT --vcf {VCF} --projectDir {project_dir}
    """.format(VCF=VCF, GENOME=GENOME, MASK=MASK, project_dir=project_dir, MU=MU)
    print(spec)
    return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)


### Function calls


for p in population_l:
    os.makedirs(path_to_output + p, exist_ok=True)
    print(p)

gwf.map(relernn, population_l, name="relernn",
                extra={"input_dir": path_to_input,
                       "output_dir": path_to_output, "GENOME": genome_sizes,
                       "MASK": mask_bed, "MU": mu})
