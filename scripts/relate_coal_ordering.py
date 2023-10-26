# Loading packages
import pandas as pd
import argparse
import tskit
import os

# Arg definition
parser = argparse.ArgumentParser(description="Analyze coal ordering from tree sequence")
parser.add_argument("-t", help="input tree")
parser.add_argument("-p", help="poplabel path")
parser.add_argument("-i", help="focal ID")
parser.add_argument("-o", help="name and path to write the finished file")


# Function for coal ordering
def coalescence_ordering(tree, IDs, sample_counts):
    df_list = []
    for i in IDs:
        pop_list = []
        gen_list = []
        current_node = i
        while tree.depth(current_node) > 0:
            # Find parent node
            parent_node = tree.parent(current_node)
            c_samples = pd.Series([x for x in tree.samples(current_node)])
            c_sample_counts = c_samples.map(i_mapping).value_counts()
            # If pop is not already added to list, add pop and note coal time when it surpasses 50 %
            for p in c_sample_counts.index:
                if p  not in pop_list and c_sample_counts[p] > sample_counts[p]/2:
                    pop_list.append(p)
                    gen_list.append(tree.time(current_node))
            current_node = parent_node
        d = {"ID": i, "sites": tree.num_sites, "span": tree.span, "start": tree.interval[0]}
        for i in range(len(pop_list)):
            d["coal_{}".format(i)] = pop_list[i]
        for i in range(len(gen_list)):
            d["coal_date_{}".format(i)] = gen_list[i]
        df_list.append(pd.DataFrame(d, index=[i]))
    return pd.concat(df_list)


# Parsing args and constant defs
args = parser.parse_args()
tree = args.t
poplabel_path = args.p
basepath = os.path.dirname(args.t)
treename = os.path.basename(args.t)


# Loading poplabels
poplabels = pd.read_csv(poplabel_path, sep=" ")

# Setup based on poplabels
if "hap" in tree:
    ID_list = list(poplabels.loc[poplabels.ID == args.i].index.astype(int))
    sample_counts = poplabels["GROUP"].value_counts()
    i_mapping = {}
    for i, row in poplabels.iterrows():
        i_mapping[i] = row.GROUP
        if row.ID == args.i:
            ID_list = [i]
else:
    sample_counts = poplabels["GROUP"].value_counts()*2
    i_mapping = {}
    for i, row in poplabels.iterrows():
        i_mapping[i*2] = row.GROUP
        i_mapping[i*2+1] = row.GROUP
        if row.ID == args.i:
            ID_list = [i*2, i*2+1]

print(ID_list)

# Loading and running through the trees
c = 0
df_list = []
ts = tskit.load(tree)

print("Starting coalescence ordering with")
print("{} ind in focus and the following sample counts {}".format(ID_list, sample_counts))
for tree in ts.trees(sample_lists=True):
    df = coalescence_ordering(tree, ID_list, sample_counts)
    if c % 2500 == 0:
        print(c)
    c += 1
    df_list.append(df)
full_df = pd.concat(df_list)

full_df.to_csv("{}/analysis_chunks/{}_{}.txt".format(basepath, args.i, treename[:-6]), index=False)
