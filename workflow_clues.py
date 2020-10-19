from gwf import Workflow, AnonymousTarget
import os
import pandas as pd
from groups import Group

gwf = Workflow()

# Absolute paths to various parts - possibly rewrite to use relative paths
path_to_relate = "/faststorage/home/eriks/relate-clues/relate_v1.1.2_x86_64_dynamic/"
path_to_prepared = "/faststorage/home/eriks/relate-clues/steps/"
path_to_script = "/faststorage/home/eriks/relate-clues/scripts/clues_master.py"
path_to_results = "/faststorage/home/eriks/relate-clues/results/"
path_to_relate_results = "/faststorage/home/eriks/relate-clues/results/"
summarize_script = "/faststorage/home/eriks/relate-clues/scripts/summarize_clues.py"
# Name of directory, followed by mutation rate and populaiton size.
pop_information = {"CEU_kept_prepared/": ["1.25e-8", "6000"],
                   "YRI_kept_prepared/": ["1.25e-8", "32000"]}
# specification of clues variables. Timebins indicates the time to test for selection
#  burn_in and thinning are used for MCMC sampling
timebins = "/faststorage/home/eriks/relate-clues/data/clues_supporting/time0_500.txt"
burn_in = 50
thinning = 10
# Chromosome number, start and end of window, as well as number of tests in total and number of tests per job.
chromosome = "3"
window_start = 45000000
window_end = 55000000
number_of_tests = 11
tests_per_job = 10
# Name added to dir
dir_suffix = "chrom{}_{}start_{}end_{}tests".format(chromosome, window_start, window_end, number_of_tests)


def relate_clues(chunk, chunk_count, script, number, out_dir, pop_inf, input_dir, relate_path, timebins, burnin, thin):
    i = input_dir+"chrom{}_popsize.coal".format(number)
    dist = i[:-5]+".dist"
    coal = input_dir+"chrom{}_popsize.coal".format(number)
    tmp = out_dir+"tmp/chrom{}_{}_branch_lengths".format(number, snp)
    o = out_dir+"tmp/chrom{}_snp{}_clues".format(number, snp)
    m = pop_inf[0]
    inputs = [i]
    outputs = o+"_output.txt"
    options = {
        "cores": 1,
        "memory": "4g",
        "walltime": "3:00:00"
    }
    spec = """
    python {} {} {} -o {} -p {} -
    """.format(script, chunk, chunk_count, )
    print(spec)
    #return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)


def summarize_clues(clues_results, script, output):
    inputs = clues_results
    i = " ".join(map(str, clues_results))
    print(i)
    outputs = output
    options = {
        "cores": 1,
        "memory": "4g",
        "walltime": "3:00:00"
    }
    spec = """
    python {} -i {} -o {}
    """.format(script, i, output)
    print(spec)
    return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)


site_spacing = (window_end-window_start)/number_of_tests
chunk_number = number_of_tests//tests_per_job + (number_of_tests % tests_per_job > 0)
chunk_list = list(range(0, job_number+1))
print("""Trying to place a site every {} bases, leading to {} tests spread across {} jobs.
      """.format(site_spacing, number_of_tests, job_number))
for pop in pop_information:
    pop_name = pop[:-10]
    out_dir = "results/{}_{}_clues/".format(pop_name, dir_suffix)
    os.makedirs(out_dir+"tmp/", exist_ok=True)
    freq = pd.read_table(path_to_relate_results+"{}_relate/chrom{}_selection.freq".format(pop_name, chromosome), sep='\s+')
    lin = pd.read_table(path_to_relate_results+"{}_relate/chrom{}_selection.lin".format(pop_name, chromosome), sep='\s+')
    sele = pd.read_table(path_to_relate_results+"{}_relate/chrom{}_selection.sele".format(pop_name, chromosome), sep='\s+')
    freq = pd.concat([freq["pos"], freq["DataFreq"], lin["0.000000"], sele["when_mutation_has_freq2"]],
                     axis=1, keys=["pos", "DataFreq", "lineages", "relate_p"])
    positions = freq.loc[(freq["pos"] >= window_start) & (freq["pos"] <= window_end)]
    next_site = window_start
    sites = []
    kept_rows = []
    for site, row in positions.iterrows():
        if row["pos"] >= next_site:
            next_site += site_spacing
            sites.append(int(row["pos"]))
    df = positions[positions.pos.isin(sites)]
    df.to_csv(out_dir+"clues_table.txt", index=False, sep=" ")
    with Group(gwf, suffix=pop_name) as g:
        total_workflow = g.map(relate_clues, chunk_list, name="chrom_{}_{}".format(chromosome, window_start),
                               extra={"chunk_count": chunk_number, "script": 
                                      "number": chromosome, "out_dir": path_to_results+"{}_clues/".format(pop_name),
                                      "pop_inf": pop_information[pop], "input_dir": path_to_prepared+"{}_relate/".format(pop_name),
                                      "relate_path": path_to_relate, "timebins": timebins,
                                      "burnin": burn_in, "thin": thinning})
        # g.target_from_template("summarize",
        #                        summarize_clues(total_workflow.outputs, script=summarize_script,
        #                                        output=path_to_results+"{}_clues/chrom{}".format(pop_name, chromosome)))

# gwf.map(clues, snp=sites, extra={"number": chromosome, "out_dir": path_to_results+"CEU_kept_clues/",
#                                  "input_dir": path_to_prepared+"CEU_kept_relate/",
#                                  "timebins": timebins})
