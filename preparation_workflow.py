from gwf import Workflow, AnonymousTarget
import os, re
import numpy as np

gwf = Workflow()

path_to_relate = "/faststorage/home/eriks/relate-clues/relate_v1.1.2_x86_64_dynamic/"
path_to_vcfs = "/faststorage/project/simons/data/1000Genomes/"
path_to_ancestor = "/faststorage/project/simons/data/1000Genomes/ancestral_alignmens/human_ancestor_GRCh37_e59/"
haps_sample_dir = "steps/haps_sample"

vcfs = []
ancestors = {}
chromosomes = list(range(1, 23))+["X"]

for chrom in chromosomes:
	if chrom == 'X':
		vcf_path_and_name = os.path.join(path_to_vcfs, 'ALL.chrX.phase3_shapeit2_mvncall_integrated_v1b.20130502.genotypes') #dont use file extension
	else:
		vcf_path_and_name = os.path.join(path_to_vcfs, 'ALL.chr{}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes'.format(chrom))
	ancestor_name = os.path.join(path_to_ancestor, "human_ancestor_{}.fa".format(chrom))
	vcfs.append({"vcf_path": vcf_path_and_name, "chrom": chrom})
	ancestors[chrom] = ancestor_name

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
	return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)

os.makedirs(haps_sample_dir, exist_ok=True)
gwf.map(vcf_to_haps, vcfs, extra = {"relate_path": path_to_relate, "output_dir": haps_sample_dir})