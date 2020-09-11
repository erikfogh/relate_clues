from gwf import Workflow, AnonymousTarget
from gwf.workflow import collect
import os, re, glob
import numpy as np
from groups import Group

gwf = Workflow()

path_to_relate = "/faststorage/home/eriks/relate-clues/relate_v1.1.2_x86_64_dynamic/"
path_to_vcfs = "/faststorage/project/simons/data/1000Genomes/"
path_to_ancestor = "/faststorage/project/simons/data/1000Genomes/ancestral_alignmens/human_ancestor_GRCh37_e59/"
haps_sample_dir = "steps/haps_sample"
pop_files = "data/pops/"
use_all = True
pop_list = glob.glob(pop_files + "*.txt")

vcfs = []
ancestors = []
chromosomes = list(range(1, 23))+["X"]

for chrom in chromosomes:
	if chrom == 'X':
		vcf_path_and_name = os.path.join(path_to_vcfs, 'ALL.chrX.phase3_shapeit2_mvncall_integrated_v1b.20130502.genotypes') #dont use file extension
	else:
		vcf_path_and_name = os.path.join(path_to_vcfs, 'ALL.chr{}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes'.format(chrom))
	ancestor_name = os.path.join(path_to_ancestor, "human_ancestor_{}.fa".format(chrom))
	vcfs.append({"vcf_path": vcf_path_and_name, "chrom": chrom})
	ancestors.append(ancestor_name)

def vcf_to_haps(vcf_path, chrom, relate_path, output_dir):
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
	print(outputs)
	return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)

def prepare_input(haps, sample, ancestor, number, pop, output_dir, relate_path):
	inputs = {"haps": haps, "sample": sample}
	PrepareInputFiles = os.path.join(relate_path, "scripts/PrepareInputFiles/PrepareInputFiles.sh")
	destination_name = output_dir + "chrom{}".format(number)
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
	{}  --haps {} --sample {} --ancestor {}{} -o {}
	""".format(PrepareInputFiles, haps, sample, ancestor, remove_ids, destination_name)
	print(outputs)
	return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)

def cleanup(files, dir):
	inoputs = [files]
	outputs = []
		options = {
		"cores": 2,
		"memory": "4g",
		"walltime": "1:00:00"
	}
	spec = """
	rm -r {}
	""".format(dir)
	return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)

os.makedirs(haps_sample_dir, exist_ok=True)
haps_sample = gwf.map(vcf_to_haps, vcfs, extra = {"relate_path": path_to_relate, "output_dir": haps_sample_dir})
l_d = haps_sample.outputs
for d in range(len(l_d)):
	l_d[d]["ancestor"] = ancestors[d]
	l_d[d]["number"] = str(chromosomes[d])


if use_all == True:
	pop_list.append("all_individuals")
prep_list = []
#print(l_d)
for pop in pop_list:
	if os.path.basename(pop).endswith("txt"):
		popname = os.path.basename(pop)[:-4]
	else:
		popname = pop
	path_to_dir = "steps/{}_prepared/".format(popname)
	os.makedirs(path_to_dir, exist_ok=True)
	with Group(gwf, suffix = popname) as g:
		preps = g.map(prepare_input, l_d, extra = {"pop": pop, "output_dir": path_to_dir, "relate_path": path_to_relate})
		prep_list.append(preps.outputs)

removal = gwf.target_from_template(
	name = "cleanup",
	cleanup(files = collect(prep_list[0], ["haps"]),
	dir = haps_sample_dir
)

print(prep_list[0])