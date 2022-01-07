import pandas as pd
import numpy as np

input_table = pd.read_csv(snakemake.input[0], delim_whitespace=True)


output_table = pd.DataFrame.from_dict({'MARKERNAME': input_table['variant.id'],
                                       'CHROMOSOME': input_table['chr'].str.slice(3),
                                       'POSITION': input_table['pos'],
                                       'N': input_table['n.obs'],
                                       'EAF': input_table['freq'],
                                       'EA': input_table['effect.allele'],
                                       'NEA': input_table['other.allele']})

if 'n.cases' in input_table.columns:
    output_table['OR'] = np.exp(input_table['Est'])
    output_table['OR_95L'] = np.exp(input_table['Est'] - 1.96*input_table['Est.SE'])
    output_table['OR_95U'] = np.exp(input_table['Est'] + 1.96 * input_table['Est.SE'])
    output_table = output_table.loc[(output_table['OR_95L'] > 0) & (output_table['OR_95U'] != np.inf)]
else:
    output_table['BETA'] = input_table['Est']
    output_table['SE'] = input_table['Est.SE']

output_table.to_csv(snakemake.output[0], sep=" ", index=False)
