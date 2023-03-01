# Loading packages
import pandas as pd
import argparse
import tskit

# Arg definition
parser = argparse.ArgumentParser(description="Analyze coal ordering from tree sequence")
parser.add_argument("-t", help="input tree")
parser.add_argument("-i", help="focal ID")
parser.add_argument("-p", help="poplabel path")
parser.add_argument("-o", help="name and path to write the finished file")

# Function for coal ordering
def coalescence_ordering(tree, IDs, sample_counts):
    df_list = []
    for i in IDs:
        pop_list = []
        gen_list = []
        coal_counts = {}
        for p in sample_counts.index:
            coal_counts[p] = 0
        current_node = i
        while tree.depth(current_node) > 0:
            # Find parent node
            parent_node = tree.parent(current_node)
            # Determine children, and then pick alternate
            # cannot find a method for this, so I use an explicit if/else
            children = tree.children(parent_node)
            if current_node == children[0]:
                alt_node = children[1]
            else:
                alt_node = children[0]
            # Determine which populations are present under the alternate node
            alt_samples = pd.Series([x for x in tree.samples(alt_node)])
            alt_sample_counts = alt_samples.map(i_mapping).value_counts()
            for p in alt_sample_counts.index:
                coal_counts[p] += alt_sample_counts[p]
            # If pop is not already added to list, add pop and note coal time
            for p in alt_sample_counts.index:
                if p  not in pop_list and coal_counts[p] > sample_counts[p]/2:
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
focal_ID = args.i
poplabel_path = args.p
out_path = args.o

# Loading data
ts = tskit.load(tree)
poplabels = pd.read_csv(poplabel_path, sep=" ",
                        names=["i", "ID", "POP", "GROUP", "SEX"], header=0)

# Setup based on poplabels
ID_t1 = poplabels.loc[poplabels.ID == focal_ID].index.values[0]
ID_list = [ID_t1, ID_t1+1]
sample_counts = poplabels["GROUP"].value_counts()*2
i_mapping = {}
for i, row in poplabels.iterrows():
    i_mapping[i*2] = row.GROUP
    i_mapping[i*2+1] = row.GROUP

# Running through the trees
c = 0
df_list = []
for tree in ts.trees():
    df = coalescence_ordering(tree, ID_list, sample_counts)
    if c % 2500 == 0:
        print(c)
    c += 1
    df_list.append(df)
full_df = pd.concat(df_list)

full_df.to_csv(out_path + focal_ID + ".txt", index=False)