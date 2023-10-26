from gwf import Workflow, AnonymousTarget
import os
from groups import Group
os.environ['NUMEXPR_MAX_THREADS'] = '8' # to suppress a warning

gwf = Workflow(defaults={"account": "baboondiversity"})

# Absolute paths to various parts - possibly rewrite to use relative paths
path_to_relate = "/home/eriks/baboondiversity/people/eriks/second_analysis_baboons/relate/relate_v1.1.7_x86_64_dynamic/"
path_to_script = "/faststorage/home/eriks/relate-clues/scripts/clues_master.py"
path_to_results = "/home/eriks/baboondiversity/people/eriks/second_analysis_baboons/results/"
path_to_relate_results = "/home/eriks/baboondiversity/people/eriks/second_analysis_baboons/steps/{}_relate/"
summarize_script = "/faststorage/home/eriks/relate-clues/scripts/summarize_clues.py"
table_gen_script = "/faststorage/home/eriks/relate-clues/scripts/table_generation.py"
# Mutation rate and population size
pop_information = ["0.57e-8", "50000"]
pop_information_x = ["0.46e-8", "25000"]
# specification of clues variables. Timebins indicates the time to test for selection
#  burn_in and thinning are used for MCMC sampling
# timebins = "/faststorage/home/eriks/relate-clues/data/clues_supporting/time0_500.txt"
relate_mcmc = 150
burn_in = 50
thinning = 1
daf_bound = "0.10_0.90"
# Chromosome number, start and end of window, as well as number of tests in total and number of tests per job.
run_combinations = [{"pop_name": "Tanzania_Olive", "chromosome": ["1", "2"],
                     "window_start": 0, "window_end": "max",
                     "number_of_tests": 100, "tests_per_job": 25,
                     "doublestack": True, "run_name": "tanz_aut_test"}]
# chromosome = "hapX"
# window_start = 0
# window_end = 144000000
# number_of_tests = 100
# tests_per_job = 10
# Name added to dir
# identifier = "prio"
# dir_suffix = "chrom{}_{}_{}_{}_{}".format(chromosome, window_start, window_end, number_of_tests, identifier)


def table_gen(run_name, script, relate_results, start, end, number_of_tests, daf_bound, doublestack, out_dir, chrom):
    i = relate_results+".freq"
    inputs = [i]
    outputs = out_dir+"{}_clues_table_start_{}.txt".format(chrom, run_name)
    options = {
        "cores": 2,
        "memory": "10g",
        "walltime": "1:00:00"
    }
    spec = """
    python {} {} {} {} {} {} {} -o {}
    """.format(script, relate_results, start, end, int(number_of_tests), daf_bound, doublestack, outputs)
    print(outputs)
    return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)


def relate_clues(chunk, chunk_number, run_name, script, number, out_dir, pop_inf, input_dir, relate_path, mcmc, burnin, inputs):
    i = input_dir+"chrom{}_popsize.coal".format(number)
    chrom = "chrom{}".format(number)
    inputs = inputs
    out_name = "snp_results/{}_chunk{}_{}_table.txt".format(chrom, chunk, run_name)
    outputs = out_dir+out_name
    walltime = 2*(tests_per_job//5)+2
    options = {
        "cores": 5,
        "memory": "40g",
        "walltime": "{}:00:00".format(walltime)
    }
    spec = """
    cd /home/eriks/faststorage/relate-clues/clues/
    python {} {} {} -i {} -o {} -n {} -m {} -r {} -t {} -c {} --mcmc {} --burnin {}
    """.format(script, chunk, chunk_number, i[:-5], out_dir, out_name, pop_inf[0], relate_path, inputs, chrom, mcmc, burnin)
    print(spec)
    return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)


def summarize_clues(clues_results, script, out_dir, table):
    inputs = clues_results
    outputs = out_dir+"clues_table.txt"
    options = {
        "cores": 2,
        "memory": "4g",
        "walltime": "1:00:00"
    }
    spec = """
    python {} -i {} -t {} -o {}
    """.format(script, ",".join(inputs), table, outputs)
    print(inputs, spec)
    return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)


# Workflow in which the targets are defined by the table_gen workflow.

for run in run_combinations:
    for chromosome in run["chromosome"]:
        pop_name = run["pop_name"]
        window_start = run["window_start"]
        window_end = run["window_end"]
        number_of_tests, tests_per_job = run["number_of_tests"], run["tests_per_job"]
        doublestack = run["doublestack"]
        run_name = run["run_name"]
        chunk_number = number_of_tests//tests_per_job + (number_of_tests % tests_per_job > 0)
        chunk_list = list(range(1, chunk_number+1))
        out_dir = os.getcwd()+"/results/{}_clues/".format(pop_name)
        if chromosome in ["hapX", "X"]:
            pop_inf = pop_information_x
        else:
            pop_inf = pop_information
        print(out_dir)
        os.makedirs(out_dir+"snp_results/", exist_ok=True)
        relate_results_path = "results/{}_relate/chrom{}_selection".format(pop_name, chromosome)
        with Group(gwf, suffix="{}".format(pop_name)) as g:
            table_clues = g.target_from_template("table_picks_{}_{}".format(chromosome, run_name),
                                                 table_gen(run_name, table_gen_script, relate_results_path, window_start,
                                                           window_end, number_of_tests, daf_bound, doublestack, out_dir,
                                                           chromosome))
            total_workflow = g.map(relate_clues, chunk_list,
                               name="chrom_{}_{}".format(chromosome, run_name),
                               extra={"run_name": run_name, "chunk_number": chunk_number, "script": path_to_script,
                                      "number": chromosome, "out_dir": out_dir, "pop_inf": pop_inf,
                                      "input_dir": path_to_relate_results.format(pop_name),
                                      "relate_path": path_to_relate,
                                      "mcmc": relate_mcmc, "burnin": burn_in,
                                      "inputs": table_clues.outputs})
            g.target_from_template("summarize_{}_{}".format(chromosome, run_name),
                               summarize_clues(total_workflow.outputs, script=summarize_script,
                                               out_dir=out_dir+"{}_{}".format(chromosome, run_name), table=table_clues.outputs))


# Alternate version with a custom table defined by hand.