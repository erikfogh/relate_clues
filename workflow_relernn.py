from gwf import Workflow, AnonymousTarget
import os
from groups import Group

gwf = Workflow()

# Input paths
# Uses the recomb environment.

path_to_input = "/home/eriks/baboondiversity/data/PG_panu3_phased_chromosomes_4_7_2021/relernn_input/"
path_to_output = "steps/relernn/"
genome_sizes = path_to_input + "autosome_sizes.bed"
mask_bed = path_to_input + "autosomes.sorted.bed"
mu = "0.57e-8"
#  chromosomes = list(range(1, 21))

population_l = ["Dendro", "Arusha", "Ngorongoro"]

# Gwf functions


def relernn(name, input_dir, output_dir, GENOME, MASK, MU):
    project_dir = output_dir + name + "/"
    output = output_dir + name + "/" + name + ".PREDICT.txt"
    VCF = input_dir + name + ".vcf"
    inputs = [VCF, GENOME, MASK]
    outputs = [output]
    options = {
        "cores": 2,
        "memory": "15g",
        "walltime": "60:00:00",
        "account": "baboondiversity"
    }
    spec = """
    ReLERNN_SIMULATE --vcf {VCF} --genome {GENOME} \
    --mask {MASK} --projectDir {project_dir} --assumedMu {MU} \
    -l 11 -r 5 --nCPU 2
    ReLERNN_TRAIN --projectDir {project_dir} --nEpochs 200
    ReLERNN_PREDICT --vcf {VCF} --projectDir {project_dir}
    """.format(VCF=VCF, GENOME=GENOME, MASK=MASK, project_dir=project_dir, MU=MU)
    return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)


def relernn_correct(name, output_dir):
    project_dir = output_dir + name + "/"
    input = output_dir + name + "/" + name + ".PREDICT.txt"
    output = output_dir + name + "/" + name + ".PREDICT.BSCORRECT.txt"
    inputs = [input]
    outputs = [output]
    options = {
        "cores": 4,
        "memory": "30g",
        "walltime": "10:00:00",
        "account": "baboondiversity"
    }
    spec = """
    ReLERNN_BSCORRECT  --projectDir {project_dir}
    """.format(project_dir=project_dir)
    print(spec)
    return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)


# Function calls


for p in population_l:
    os.makedirs(path_to_output + p, exist_ok=True)
    print(p)

gwf.map(relernn, population_l, name="relernn",
        extra={"input_dir": path_to_input,
               "output_dir": path_to_output, "GENOME": genome_sizes,
               "MASK": mask_bed, "MU": mu})

gwf.map(relernn_correct, population_l, name="relernn_correct",
        extra={"output_dir": path_to_output})
