from gwf import Workflow, AnonymousTarget
from gwf.workflow import collect
import os
import glob
from groups import Group

gwf = Workflow()

# Associated variables and setup which can be added to environment variable instead.
path_to_relate = "/faststorage/home/eriks/relate-clues/relate_v1.1.2_x86_64_dynamic/"
path_to_vcfs = "/faststorage/project/simons/data/1000Genomes/"
path_to_ancestor = "/faststorage/project/simons/data/1000Genomes/ancestral_alignmens/human_ancestor_GRCh37_e59/human_ancestor_{}.fa"
path_to_poplabels = "data/pops/1000GP_Phase3.sample"
path_to_mask = "data/depth_chr_masks/20140520.chr{}.strict_mask.fasta"
maskx = "data/depth_chr_masks/20141020.chrX.strict_mask.fasta"
genetic_map = "/faststorage/home/eriks/relate-clues/data/recombination_maps/genetic_map_chr{}_combined_b37.txt"
genetic_map_x = "/faststorage/home/eriks/relate-clues/data/recombination_maps/genetic_map_chrX_nonPAR_combined_b37.txt"
haps_sample_dir = "steps/haps_sample"

pop_files = "data/pops/"
use_all = True
pop_list = glob.glob(pop_files + "*.txt")

if use_all is True:
    pop_list.append("all_individuals")

vcfs = []
chromosomes = list(range(1, 23))+["X"]

for chrom in chromosomes:
    if chrom == 'X':
        vcf_path_and_name = os.path.join(
            path_to_vcfs, 'ALL.chrX.phase3_shapeit2_mvncall_integrated_v1b.20130502.genotypes'
            )  # dont use file extension
    else:
        vcf_path_and_name = os.path.join(
            path_to_vcfs, 'ALL.chr{}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes'.format(chrom)
            )
    ancestor_name = os.path.join(path_to_ancestor, "human_ancestor_{}.fa".format(chrom))
    vcfs.append({"vcf_path": vcf_path_and_name, "chrom": chrom})


def vcf_to_haps(vcf_path, chrom, relate_path, output_dir):
    """Converts vcf files to haps/sample using the script from Relate """
    inputs = [vcf_path+".vcf.gz"]
    haps_out = os.path.join(output_dir, "chrom{}.haps".format(chrom))
    sample_out = os.path.join(output_dir, "chrom{}.sample".format(chrom))
    RelateFileFormats = os.path.join(relate_path, "bin/RelateFileFormats")
    outputs = {"haps": haps_out, "sample": sample_out}
    options = {
        "cores": 2,
        "memory": "4g",
        "walltime": "1:00:00"
    }
    spec = """
    {} --mode ConvertFromVcf --haps {} --sample {} -i {}
    """.format(RelateFileFormats, haps_out, sample_out, vcf_path)
    return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)


def prepare_input(haps, sample, mask, ancestor, pop, output_dir, poplabels, relate_path):
    inputs = {"haps": haps, "sample": sample}
    s = os.path.basename(haps)
    number = s[5:s.find(".")]
    PrepareInputFiles = os.path.join(relate_path, "scripts/PrepareInputFiles/PrepareInputFiles.sh")
    destination_name = output_dir + "chrom{}".format(number)
    if number == "X":
        n_mask = maskx
    else:
        n_mask = mask.format(number)
    ancestor_in = ancestor.format(number)
    outputs = {"haps": destination_name+".haps.gz", "sample": destination_name+".sample.gz"}
    if pop != "all_individuals":
        remove_ids = " --remove_ids " + pop
    else:
        remove_ids = ""
    options = {
        "cores": 8,
        "memory": "4g",
        "walltime": "4:00:00"
    }
    spec = """
    {}  --haps {} --sample {} --mask {} --ancestor {}{} --poplabels {} -o {}
    """.format(PrepareInputFiles, haps, sample, n_mask, ancestor_in, remove_ids, poplabels, destination_name)
    return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)


def make_chunks(haps, sample, pop, output_dir, relate_path):
    inputs = {"haps": haps, "sample": sample}
    s = os.path.basename(haps)
    number = s[5:s.find(".")]
    make_chunks = os.path.join(relate_path, "bin/Relate")
    dist = output_dir+"chrom{}.dist.gz".format(number)
    if number == "X":
        gene_map = genetic_map_x
    else:
        gene_map = genetic_map.format(number)
    outputs = [output_dir+"temp/chrom{}/chunk_0.hap".format(number)]
    options = {
        "cores": 4,
        "memory": "16g",
        "walltime": "4:00:00"
    }
    spec = """
    cd {}
    {}  --mode "MakeChunks" --haps {} --sample {} --map {} --dist {} -o {} --memory 14
    """.format(output_dir+"temp/",
               make_chunks, haps, sample, gene_map, dist, "chrom{}".format(number))
    print(spec)
    return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)


os.makedirs(haps_sample_dir, exist_ok=True)
haps_sample = gwf.map(vcf_to_haps, vcfs, extra={"relate_path": path_to_relate, "output_dir": haps_sample_dir})
prep_list = []


for pop in pop_list:
    if os.path.basename(pop).endswith("txt"):
        popname = os.path.basename(pop)[:-4]
    else:
        popname = pop
    path_to_dir = os.getcwd()+"/steps/{}_prepared/".format(popname)
    os.makedirs(path_to_dir, exist_ok=True)
    os.makedirs(path_to_dir+"temp/", exist_ok=True)
    with Group(gwf, suffix=popname) as g:
        preps = g.map(prepare_input, haps_sample.outputs, name="prepare_input", extra={
                        "mask": path_to_mask, "ancestor": path_to_ancestor, "pop": pop,
                        "output_dir": path_to_dir, "poplabels": path_to_poplabels,
                        "relate_path": path_to_relate})
        chunks = g.map(make_chunks, preps.outputs, name="make_chunks", extra={
                        "pop": pop,
                        "output_dir": path_to_dir,
                        "relate_path": path_to_relate})

# Try at cleaning up the haps_sample dir, does not work currently
# gwf.target_from_template(
#     name="remove_intermediate_files",
#     template=cleanup(haps=prep_list, directory=haps_sample_dir)
# )
