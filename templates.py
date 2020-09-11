from gwf import Workflow, AnonymousTarget
import os, re
import numpy as np

path_to_relate = "/faststorage/home/eriks/relate-clues/relate_v1.1.2_x86_64_dynamic/bin/"
path_to_vcfs = "/faststorage/project/simons/data/1000Genomes/"
vcf_name = "ALL.chr22.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz"

def vcf_to_haps(relate_path, haps_out, sample_out, input):
	"""RelateFileFormats to convert from vcf to haps/sample"""
	inputs = {"input": input}
	out_file = os.path.join(outpath, "chr_"+number)
	outputs = [out_file+".frq.count"]
	options = {
		"cores": 4,
		"memory": "4g"
	}
	spec = """
	vcftools --gzvcf {} --counts --min-alleles 2 --max-alleles 2 --out {} --keep {}
	""".format(path, out_file, pop)
	return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)