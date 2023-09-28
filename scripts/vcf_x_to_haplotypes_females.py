import sys
import pandas as pd

# the metadata is used to seperate samples/species.

# Input file here is a phased file only containing females.

c = 0

with sys.stdin as f:
    for line in f:
        if line.startswith('#'):
            # header
            if line.startswith('#CHROM'):
                CHROM, POS, ID, REF, ALT, QUAL, FILTER, INFO, FORMAT, *all_samples = line.split()
                female_ids = []
                for i in range(len(all_samples)):
                    female_ids.append(all_samples[i]+"_1")
                    female_ids.append(all_samples[i]+"_2")
                start_line = [CHROM, POS, ID, REF, ALT, QUAL, FILTER, INFO, FORMAT]
                start_line.extend(female_ids)
                print('\t'.join(start_line))
            else:
                print(line, end='')
        else:
            # calls
            CHROM, POS, ID, REF, ALT, QUAL, FILTER, INFO, FORMAT, *calls = line.split()
            #
            # for skipping POS in the PAR.
            if int(POS) < 2500000:
                continue
            # for only taking small chunks (testing)
            # c += 1
            # if c > 10:
            #     break
            female_calls = []
            for i in range(len(all_samples)):
                female_calls.extend([calls[i][0], calls[i][2]])
            start_line = [CHROM, POS, ID, REF, ALT, QUAL, FILTER, INFO, FORMAT]
            start_line.extend(female_calls)
            print('\t'.join(start_line))
