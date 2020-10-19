import os
import numpy as np
from gwf import Workflow, AnonymousTarget

gwf = Workflow()

path_to_vcfs = "/faststorage/project/simons/data/1000Genomes/"
outpath = "/faststorage/home/eriks/relate-clues/data/snp_lists/"

vcfs = []
chromosomes = list(range(1, 23))+["X"]

for chrom in chromosomes:
    if chrom == 'X':
        vcf_path_and_name = os.path.join(
            path_to_vcfs, 'ALL.chrX.phase3_shapeit2_mvncall_integrated_v1b.20130502.genotypes.vcf.gz'
            )  # dont use file extension
    else:
        vcf_path_and_name = os.path.join(
            path_to_vcfs, 'ALL.chr{}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz'.format(chrom)
            )
    vcfs.append({"vcf_path": vcf_path_and_name, "chrom": str(chrom)})


def snp_list(vcf_path, chrom, outpath):
    """VCFtools returning a file containing counts."""
    inputs = {"path": vcf_path}
    out_file = os.path.join(outpath, "chrom"+chrom)
    outputs = [out_file+".frq.count"]
    options = {
        "cores": 10,
        "memory": "4g",
        "walltime": "2:00:00"
    }
    spec = """
    vcftools --gzvcf {} --counts --min-alleles 2 --max-alleles 2 --out {} --max-missing 1
    """.format(vcf_path, out_file)
    return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)


os.makedirs(outpath, exist_ok=True)
gwf.map(snp_list, vcfs, extra={"outpath": outpath})
