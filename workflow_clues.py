from gwf import Workflow, AnonymousTarget
import os
from groups import Group
os.environ['NUMEXPR_MAX_THREADS'] = '8' # to suppress a warning

gwf = Workflow()

# Absolute paths to various parts - possibly rewrite to use relative paths
path_to_relate = "/faststorage/home/eriks/relate-clues/relate_v1.1.2_x86_64_dynamic/"
path_to_prepared = "/faststorage/home/eriks/relate-clues/steps/"
path_to_script = "/faststorage/home/eriks/relate-clues/scripts/clues_master.py"
path_to_results = "/faststorage/home/eriks/relate-clues/results/"
path_to_relate_results = "/faststorage/home/eriks/relate-clues/results/"
summarize_script = "/faststorage/home/eriks/relate-clues/scripts/summarize_clues.py"
table_gen_script = "/faststorage/home/eriks/relate-clues/scripts/table_generation.py"
# Name of directory, followed by mutation rate and populaiton size.
pop_information = {"CEU_kept_prepared/": ["1.25e-8", "30000"],
                   } # "YRI_kept_prepared/": ["1.25e-8", "32000"]
# specification of clues variables. Timebins indicates the time to test for selection
#  burn_in and thinning are used for MCMC sampling
timebins = "/faststorage/home/eriks/relate-clues/data/clues_supporting/time1500_2500.txt"
relate_mcmc = 150
burn_in = 50
thinning = 1
daf_bound = 0.25
prioritize_sites = False
# Chromosome number, start and end of window, as well as number of tests in total and number of tests per job.
chromosome = "2"
window_start = 136000000
window_end = 137000000
number_of_tests = 100
tests_per_job = 10
# Name added to dir
identifier = "no_prio_ooa"
dir_suffix = "chrom{}_{}_{}_{}_{}".format(chromosome, window_start, window_end, number_of_tests, identifier)


def table_gen(script, relate_results, start, end, spacing, daf_bound, prioritize, out_dir):
    i = relate_results+".freq"
    inputs = [i]
    outputs = out_dir+"clues_table_temp.txt"
    options = {
        "cores": 2,
        "memory": "10g",
        "walltime": "1:00:00"
    }
    spec = """
    python {} {} {} {} {} {} {} -o {}
    """.format(script, relate_results, start, end, int(spacing), daf_bound, prioritize, out_dir)
    print(spec)
    return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)


def relate_clues(chunk, chunk_number, script, number, out_dir, pop_inf, input_dir, relate_path, timebins, mcmc, burnin, thin):
    i = input_dir+"chrom{}_popsize.coal".format(number)
    inputs = [out_dir+"clues_table_temp.txt"]
    outputs = out_dir+"chunk{}_table.txt".format(chunk)
    walltime = 10*(tests_per_job//10)+10
    options = {
        "cores": 10,
        "memory": "30g",
        "walltime": "{}:00:00".format(walltime)
    }
    spec = """
    cd clues/
    python {} {} {} -i {} -o {} -m {} -r {} -b {} --mcmc {} --burnin {} --thin {}
    """.format(script, chunk, chunk_number, i[:-5], out_dir, pop_inf[0], relate_path, timebins, mcmc, burnin, thin)
    print(spec)
    return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)


def summarize_clues(clues_results, script, out_dir):
    inputs = clues_results
    outputs = out_dir+"clues_table.txt"
    options = {
        "cores": 2,
        "memory": "4g",
        "walltime": "1:00:00"
    }
    spec = """
    python {} -i {} -o {}
    """.format(script, out_dir, outputs)
    print(spec)
    return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)


site_spacing = (window_end-window_start)/number_of_tests
chunk_number = number_of_tests//tests_per_job + (number_of_tests % tests_per_job > 0)
chunk_list = list(range(1, chunk_number+1))
print("""Trying to place a site every {} bases, leading to {} tests.
      """.format(site_spacing, number_of_tests))
for pop in pop_information:
    pop_name = pop[:-10]
    out_dir = os.getcwd()+"/results/{}_{}_clues/".format(pop_name, dir_suffix)
    os.makedirs(out_dir+"tmp/", exist_ok=True)
    relate_results_path = path_to_relate_results+"{}_relate/chrom{}_selection".format(pop_name, chromosome)
    with Group(gwf, suffix=pop_name+"_"+identifier) as g:
        g.target_from_template("temp_table",
                               table_gen(table_gen_script, relate_results_path, window_start,
                                         window_end, site_spacing, daf_bound, prioritize_sites, out_dir))
        total_workflow = g.map(relate_clues, chunk_list,
                               name="chrom_{}_{}_{}".format(chromosome, window_start, number_of_tests),
                               extra={"chunk_number": chunk_number, "script": path_to_script,
                                      "number": chromosome, "out_dir": out_dir,
                                      "pop_inf": pop_information[pop], "input_dir": path_to_prepared+"{}_relate/".format(pop_name),
                                      "relate_path": path_to_relate, "timebins": timebins,
                                      "mcmc": relate_mcmc, "burnin": burn_in, "thin": thinning})
        g.target_from_template("summarize",
                               summarize_clues(total_workflow.outputs, script=summarize_script,
                                               out_dir=out_dir))

