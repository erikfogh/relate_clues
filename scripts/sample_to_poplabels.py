import pandas as pd

path_to_file = "data/pops/1000GP_Phase3.sample"

sample_file = pd.read_table(path_to_file, sep = " ")

poplabels_file = pd.DataFrame()

sample = sample_file["ID"]
population = sample_file["pop"]
group = sample_file["group"]
sex = []

print(sample_file.columns)
