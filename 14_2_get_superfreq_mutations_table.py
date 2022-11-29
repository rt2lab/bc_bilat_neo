import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import scipy as sp
from matplotlib_venn import venn2
import pyvenn.venn as venn
from matplotlib.patches import Patch
import collections
from matplotlib.colors import ListedColormap
from sklearn import preprocessing
from scipy.spatial.distance import pdist, squareform
import matplotlib.gridspec as gridspec
import matplotlib as mpl
import numpy as np

pd.options.display.max_columns=300
pd.options.display.max_rows=1000

import warnings
warnings.filterwarnings("ignore")


patient_list = ['P{}'.format(i) for i in range(1, 7)]
patient_side_list = ['{}_L'.format(p) for p in patient_list] + \
    ['{}_R'.format(p) for p in patient_list]
table_list = list()
purity_dict = dict()
for ps in patient_side_list:
    df = pd.read_csv('results/superFreq/plots/{}/somaticVariants.csv'.format(ps), sep=',', )
    df.columns = ['p_mutation_id'] + df.columns[1:].to_list()
    df = df.assign(mutation_id=df.chr.astype(str) + '_' + df.start.astype(int).astype(str))
    df = df.assign(mutation_type='SNV')
    df.loc[(df.reference.str.len())!=(df.variant.str.len()), 'mutation_type'] = 'indel'
    df.loc[df.start!=df.end, 'mutation_type'] = 'indel'
    river = pd.read_csv('results/superFreq/plots/{}/rivers/{}-river.tsv'.format(ps, ps), sep='\t')

    river_snv_indel = river[~river.name.str.contains('[0-9]+[kM]?bp [AB]+')].dropna(subset=['chr', 'start'], how='any')
    river_snv_indel = river_snv_indel.assign(mutation_id=river_snv_indel.chr.astype(str) + '_' + river_snv_indel.start.astype(int).astype(str))
    cols_river_id = ['mutation_id', 'clone']
    cols_river_var = [c for c in river_snv_indel.columns if 'clonality' in c]
    river_melt = pd.melt(river_snv_indel[cols_river_id + cols_river_var],
                         id_vars=cols_river_id, value_vars=cols_river_var,
                         value_name='ccf', var_name='pre_sample')
    river_melt = river_melt.assign(sample=river_melt.pre_sample.str.split('.').str[1])
    df_clone = pd.merge(df, river_melt, on=['sample', 'mutation_id'], how='left')
    table_list.append(df_clone)
    
    clonal_clone = river_snv_indel.clone.min()
    pur = river_snv_indel[river_snv_indel.clone==clonal_clone].mean(axis=0)
    purity_dict.update(pur[pur.index.str.contains('clonality')].to_dict())

purity_df = pd.DataFrame.from_dict(purity_dict, orient='index', columns=['purity'])
purity_df = purity_df.assign(sample=purity_df.index.str.replace('clonality.', ''))
purity_df = purity_df[~purity_df['sample'].str.contains('GL')]
purity_df.loc[purity_df['sample']=='RD6_R', 'sample'] = 'RD6B_R'
purity_df = purity_df.assign(patient_id=purity_df['sample'].str[2])
purity_df = purity_df.assign(side_id=purity_df['sample'].str.extract('([AB]?_[LR])'))

purity_df = purity_df.assign(side_id=pd.Categorical(purity_df.side_id, ordered=True, categories=['_L', 'A_L', 'B_L', '_R', 'A_R', 'B_R']))
purity_df = purity_df.assign(sampletype=purity_df['sample'].str[:2])

purity_df.side_id.unique()

all_mut_raw = pd.concat(table_list, axis=0)
all_mut_raw_somatic = all_mut_raw[(all_mut_raw.somaticScore>0.5)&(~all_mut_raw['sample'].str.contains('GL'))&(all_mut_raw.clone!='germline')]

all_mut_raw_somatic.loc[all_mut_raw_somatic['sample']=='RD6_R', 'sample'] = 'RD6B_R'
all_mut_raw_somatic = all_mut_raw_somatic.assign(patient_id=all_mut_raw_somatic['sample'].str[2])
all_mut_raw_somatic = all_mut_raw_somatic.assign(side_id=all_mut_raw_somatic['sample'].str.extract('([AB]?_[LR])'))
all_mut_raw_somatic = all_mut_raw_somatic.assign(sampletype=all_mut_raw_somatic['sample'].str[:2])
all_mut_raw_somatic = all_mut_raw_somatic.assign(side_id=pd.Categorical(all_mut_raw_somatic.side_id, ordered=True, categories=['_L', 'A_L', 'B_L', '_R', 'A_R', 'B_R']))
type_mut_dict = {k: i+1 for i, k in enumerate(all_mut_raw_somatic.groupby('type').severity.median().sort_values().index.to_list())}
all_mut_raw_somatic = all_mut_raw_somatic.assign(type_num=all_mut_raw_somatic['type'].map(type_mut_dict))


# get list of mutations to manually check
# those mutations were removed, and then superFreq was re-run
# with the final results, this table is empty
all_mut_raw_somatic_tweak = all_mut_raw_somatic.copy()
all_mut_raw_somatic_tweak = all_mut_raw_somatic_tweak.assign(new_sample=all_mut_raw_somatic_tweak['sample'])
all_mut_raw_somatic_tweak.loc[(all_mut_raw_somatic_tweak['sample'].str.contains('6'))&(all_mut_raw_somatic_tweak['sample'].str.contains('R')), 'new_sample'] = 'PT6_R'
all_mut_raw_somatic_tweak.loc[(all_mut_raw_somatic_tweak['sample'].str.contains('6'))&(all_mut_raw_somatic_tweak['sample'].str.contains('L')), 'new_sample'] = 'PT6_L'
print(all_mut_raw_somatic_tweak['sample'].unique())
table_list = list()
for patient in range(1, 7):
    sample_list = ['PT{}_{}'.format(patient, side) for side in ('L', 'R')]
    sub_mut_df = all_mut_raw_somatic_tweak[all_mut_raw_somatic_tweak['new_sample'].isin(sample_list)]
    sub_mut_df.drop_duplicates(subset=['new_sample', 'mutation_id'], inplace=True)
    mutc = sub_mut_df.mutation_id.value_counts()
    ii = mutc[mutc>=2].index.to_list()
    print(sub_mut_df[sub_mut_df.mutation_id.isin(ii)])
    table_list.append(sub_mut_df[sub_mut_df.mutation_id.isin(ii)])


all_mut_raw_somatic.to_csv('bilat_mutations.csv', sep='\t', index=False)




