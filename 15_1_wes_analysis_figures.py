%matplotlib inline
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import scipy as sp
from matplotlib_venn import venn2
import pyvenn.venn as venn
from matplotlib.patches import Patch
from matplotlib.colors import ListedColormap
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

all_mut_raw_somatic_filter = pd.read_csv('bilat_mutations.csv',
                                         sep='\t', index=False)

###########################
# get relevant statistics #
###########################

print("median number of mutations in PT samples", all_mut_raw_somatic_filter[all_mut_raw_somatic_filter.sampletype=='PT']['sample'].value_counts().median())
print("median number of mutations in RD samples", all_mut_raw_somatic_filter[all_mut_raw_somatic_filter.sampletype=='RD']['sample'].value_counts().median())
print("median number of mutations in all samples", all_mut_raw_somatic_filter['sample'].value_counts().median())
print("median number of driver mutations in all samples", all_mut_raw_somatic_filter[(all_mut_raw_somatic_filter.isCosmicCensus)&(all_mut_raw_somatic_filter.severity<=10)]['sample'].value_counts().median())
print("median number of driver mutations in RD samples", all_mut_raw_somatic_filter[(all_mut_raw_somatic_filter.sampletype=='RD')&(all_mut_raw_somatic_filter.isCosmicCensus)&(all_mut_raw_somatic_filter.severity<=10)]['sample'].value_counts().median())
print("median number of driver mutations in PT samples", all_mut_raw_somatic_filter[(all_mut_raw_somatic_filter.sampletype=='PT')&(all_mut_raw_somatic_filter.isCosmicCensus)&(all_mut_raw_somatic_filter.severity<=10)]['sample'].value_counts().median())

PCR_PT = ['PT1_R', 'PT2_R', 'PT2_L', 'PT4_L', 'PT5_L', 'PT6A_R', 'PT6A_L', 'PT6B_L']
noPCR_PT = ['PT1_L', 'PT3_R', 'PT3_L', 'PT4_R', 'PT5_R', 'PT6B_R']

print('median nb de mutations PT PCR', all_mut_raw_somatic_filter[all_mut_raw_somatic_filter['sample'].isin(PCR_PT)]['sample'].value_counts().median())
print('median nb de mutations PT no PCR', all_mut_raw_somatic_filter[all_mut_raw_somatic_filter['sample'].isin(noPCR_PT)]['sample'].value_counts().median())
print('median nb de mutations "HFI" PT PCR', all_mut_raw_somatic_filter[(all_mut_raw_somatic_filter['sample'].isin(PCR_PT))&(all_mut_raw_somatic_filter.isCosmicCensus)&(all_mut_raw_somatic_filter.severity<=10)]['sample'].value_counts().median())
print('median nb de mutations "HFI" PT no PCR', all_mut_raw_somatic_filter[(all_mut_raw_somatic_filter['sample'].isin(noPCR_PT))&(all_mut_raw_somatic_filter.isCosmicCensus)&(all_mut_raw_somatic_filter.severity<=10)]['sample'].value_counts().median())

TNBC_PT = ['PT1_R', 'PT1_L', 'PT2_R', 'PT2_L', 'PT5_R']
luminal_PT = ['PT3_R', 'PT3_L', 'PT4_R', 'PT4_L', 'PT5_L', 'PT6A_R', 'PT6B_R', 'PT6A_L', 'PT6B_L']

print('median nb de mutations PT TNBC', all_mut_raw_somatic_filter[all_mut_raw_somatic_filter['sample'].isin(TNBC_PT)]['sample'].value_counts().median())
print('median nb de mutations PT luminal', all_mut_raw_somatic_filter[all_mut_raw_somatic_filter['sample'].isin(luminal_PT)]['sample'].value_counts().median())
print('median nb de mutations "HFI" PT TNBC', all_mut_raw_somatic_filter[(all_mut_raw_somatic_filter['sample'].isin(TNBC_PT))&(all_mut_raw_somatic_filter.isCosmicCensus)&(all_mut_raw_somatic_filter.severity<=10)]['sample'].value_counts().median())
print('median nb de mutations "HFI" PT luminal', all_mut_raw_somatic_filter[(all_mut_raw_somatic_filter['sample'].isin(luminal_PT))&(all_mut_raw_somatic_filter.isCosmicCensus)&(all_mut_raw_somatic_filter.severity<=10)]['sample'].value_counts().median())


ratio_list = list()
for ii, side in enumerate(['1_L', '3_L', '3_R', '5_R']):
    sample_list = ['{}{}'.format(sampletype, side) for sampletype in ('PT', 'RD')]
    sub_mut_df = all_mut_raw_somatic_filter[all_mut_raw_somatic_filter['sample'].isin(sample_list)]
    PT_mut = all_mut_raw_somatic_filter[all_mut_raw_somatic_filter['sample']=='PT{}'.format(side)].mutation_id.to_list()
    RD_mut = all_mut_raw_somatic_filter[all_mut_raw_somatic_filter['sample']=='RD{}'.format(side)].mutation_id.to_list()
    union_pt_rd = len(set(PT_mut).union(set(RD_mut)))
    intersection_pt_rd = len(set(PT_mut).intersection(set(RD_mut)))
    ratio_list.append( 2 * intersection_pt_rd / union_pt_rd * 100)
print("median percentage of shared mutations between PT and RD", np.median(ratio_list))

all_mut_raw_somatic_filter[all_mut_raw_somatic_filter.inGene=='TP53']

all_mut_raw_somatic_filter[all_mut_raw_somatic_filter.inGene=='CDH1']

all_mut_raw_somatic_filter[all_mut_raw_somatic_filter.inGene=='EOMES']

all_mut_raw_somatic_filter['sample'].value_counts()

####################
# figure 2 Panel A #
####################
tri_context_mut = pd.read_csv('bilat_mutations_tricontext.csv', sep='\t')
tri_context_mut = tri_context_mut.assign(nt_change=tri_context_mut.apply(lambda x: str(x['tricontext'])[2:5], axis=1))
tri_context_mut.loc[tri_context_mut.nt_change=='n', 'nt_change'] = 'indel'
uu = pd.crosstab(tri_context_mut.nt_change, tri_context_mut['sample'], normalize='columns').T
clin_data = pd.read_csv('data_clin_min.csv', sep='\t')
clin_data.index = clin_data.samplename

sns.set_style('white')
sns.set_context('poster', font_scale=1)

cosmic_only = all_mut_raw_somatic_filter[((all_mut_raw_somatic_filter.isCosmicCensus) | (
    all_mut_raw_somatic_filter.ClinVar_ClinicalSignificance.str.contains('ath'))) & (all_mut_raw_somatic_filter.severity <= 10)]
piv_order = pd.pivot_table(cosmic_only, index='patient_id',
                     columns='inGene', values='type_num', aggfunc='min')
piv = pd.pivot_table(cosmic_only, index='sample',
                     columns='inGene', values='type_num', aggfunc='min')
order_cols = (piv_order > 0).sum(axis=0).sort_values(ascending=False).index

order_row = ["PT1_R", "PT1_L", "RD1_L",
             "PT2_R", "PT2_L",
             "PT3_R", "RD3_R", "PT3_L", "RD3_L",
             "PT4_R", "PT4_L",
             "PT5_R", "RD5_R", "PT5_L",
             "PT6A_R", "PT6B_R", "PT6A_L", "PT6B_L"]



fig = plt.figure(figsize=(40, 10))
gs = gridspec.GridSpec(1, 4, width_ratios=[0.1, 0.6, 0.15, 0.15])
gs.update(wspace=0.01, hspace=0.01)
ax_main = plt.subplot(gs[0, 1])
ax_mut_count = plt.subplot(gs[0, 3])  # , sharey=ax_main)
ax_mut_type = plt.subplot(gs[0, 2])  # , sharey=ax_main)
ax_heat = plt.subplot(gs[0, 0])

cmap = sns.color_palette("Paired", 4) 
sns.heatmap(piv.loc[order_row][order_cols], linewidths=.5,
            linecolor='gray', ax=ax_main, cmap=cmap, cbar=False)
#ax_main.xaxis.tick_top()  # x axis on top
#ax_main.xaxis.set_label_position('top')

sub_type_mut_dict = {k:type_mut_dict[k] for k in type_mut_dict if type_mut_dict[k] in (1, 2, 3, 4)}
color_type_dict = {k:cmap[sub_type_mut_dict[k] - 1] for k in sub_type_mut_dict}

ax_main.tick_params(length=0, pad=10)
ax_main.set_xticklabels(ax_main.get_xticklabels(), rotation=90, ha="center")
ax_main.set_xlabel('')
ax_main.set_ylabel('')
bottom, top = ax_main.get_ylim()
ax_main.set_ylim(bottom + 0.05, top - 0.05)

sns.countplot(y='sample', data=all_mut_raw_somatic_filter,
              ax=ax_mut_count, color='#6b6b6b', order=order_row, alpha=0.7,
             orient='h')
#ax_mut_count.get_xaxis().set_visible(False)
ax_mut_count.set_xlabel('number of mutations', rotation=0, labelpad=15)
#ax_mut_count.invert_yaxis()
#ax_mut_count.set_xticklabels( [''] + ax_mut_count.get_xticklabels()[1:])
plt.setp(ax_mut_count.get_xticklabels()[0], visible=False)

uu.loc[order_row[::-1]].plot.barh(
    stacked=True, ax=ax_mut_type, alpha=0.7, width=0.7, cmap='RdYlBu')
ax_mut_type.invert_xaxis()
ax_mut_type.get_yaxis().set_visible(False)
# ax_mut_type.get_legend().remove()
ax_mut_type.legend(bbox_to_anchor=(2, 0.6), loc=2, frameon=False, title="substitution type")
ax_mut_type.set_xlabel('fraction of mutations', rotation=0, labelpad=15)



TNBC_PT = ['PT1_R', 'PT1_L', 'PT2_R', 'PT2_L', 'PT5_R']
luminal_PT = ['PT3_R', 'PT3_L', 'PT4_R', 'PT4_L',
              'PT5_L', 'PT6A_R', 'PT6B_R', 'PT6A_L', 'PT6B_L']

collevelBRCA = {"BRCA2 mut": "#0B3954",
                "BRCA1 mut": "#D64550",
                "wt": "#ADD7F6"}

colcouple_PT_RD = {"PT_with_pCR": "#49AAE3",
                   "PT_with_RD": "#F5AE78",
                   "RD": "#57EBC3"}

collevelsubtype = {"TNBC": "#D55E00",
                   "luminal": "#0072B2"}


coltumortype = {"PT": "#49AAE3",
                "RD": "#57EBC3"}

collevelside = {"R": "#D6F8D6",
                "L": "#9DC3C2"}

allcollevelpatient  = {"Patient1": "#FCB711",
                        "Patient2": "#0089D0",
                        "Patient3": "#CC004C",
                        "Patient4": "#6460AA",
                        "Patient5": "#F37021",
                        "Patient6": "#0DB14B"}

clin_data = clin_data.assign(
    color_subtype=clin_data.subtype.map(collevelsubtype))
clin_data = clin_data.assign(
    color_tumortype2=clin_data.tumortype2.map(colcouple_PT_RD))
clin_data = clin_data.assign(color_brca=clin_data.BRCA1_2_mut.map(collevelBRCA))
clin_data = clin_data.assign(color_type=clin_data.tumortype.map(coltumortype))
clin_data = clin_data.assign(color_side=clin_data.side.map(collevelside))
clin_data = clin_data.assign(color_patient=clin_data.patient.map(allcollevelpatient))

colors = [clin_data.loc[order_row].color_subtype.values,
          clin_data.loc[order_row].color_brca.values,
          clin_data.loc[order_row].color_side.values,
          clin_data.loc[order_row].color_type.values,
          clin_data.loc[order_row].color_patient.values]

m, n = 5, len(order_row)
matrix = np.zeros((m, n), int)
unique_colors = {}
for i, inner in enumerate(colors):
    for j, color in enumerate(inner):
        idx = unique_colors.setdefault(color, len(unique_colors))
        matrix[i, j] = idx
cmap2 = mpl.colors.ListedColormap(list(unique_colors))
sns.heatmap(matrix.T, cmap=cmap2, cbar=False, ax=ax_heat,
            yticklabels=order_row, xticklabels=['subtype', 'BRCA status', 'side', 'sample type', 'patient'])
#ax_heat.xaxis.tick_top() # x axis on top
#ax_heat.xaxis.set_label_position('top')
ax_heat.tick_params(length=0)
ax_heat.set_xticklabels(ax_heat.get_xticklabels(), rotation=90, ha="center")


dict_list = {'subtype': collevelsubtype,
             'BRCA status': collevelBRCA,
             'side': collevelside,
             'sample type': coltumortype,
             'patient': allcollevelpatient}

patchList = []
for dict_name, color_dict in dict_list.items():
    patchList.append(Patch(color='white', label=dict_name + ''))
    for key in color_dict:
        data_key = Patch(color=color_dict[key], label=key)
        patchList.append(data_key)
ax_heat.legend(handles=patchList, loc=2, bbox_to_anchor=(-1.5, 1.05), frameon=False )

patchList2 = []
patchList2.append(Patch(color='white', label='mutation type'))
for key in color_type_dict:
    data_key = Patch(color=color_type_dict[key], label=key)
    patchList2.append(data_key)
ax_mut_count.legend(handles=patchList2, loc=2, bbox_to_anchor=(1, 1),
                    frameon=False, title='mutation type')

sep_list = [3, 5, 9, 11, 14]
for y in sep_list:
    ax_main.axhline(y=y, color='black', lw=5)
    ax_heat.axhline(y=y, color='black', lw=5)
    ax_mut_count.axhline(y=y-0.5, color='black', lw=5)
for y in [4, 7, 9, 13, 15]:
    ax_mut_type.axhline(y=y-0.5, color='black', lw=5)
plt.savefig('figures/fig2_panelA.pdf', bbox_inches="tight")





####################
# figure 2 Panel B #
####################



fig = plt.figure(figsize=(20, 8))
gs = gridspec.GridSpec(4, 12, height_ratios=[1, 1, 2, 1])
gs.update(wspace=0.01, hspace=0.0001)
for patient in range(1, 7):
    if patient <= 5:
        sample_list = ['PT{}_{}'.format(patient, side) for side in ('R', 'L')]
        sub_mut_df = all_mut_raw_somatic_filter[all_mut_raw_somatic_filter['sample'].isin(sample_list)]
        L_mut = all_mut_raw_somatic_filter[all_mut_raw_somatic_filter['sample']=='PT{}_L'.format(patient)].mutation_id.to_list()
        R_mut = all_mut_raw_somatic_filter[all_mut_raw_somatic_filter['sample']=='PT{}_R'.format(patient)].mutation_id.to_list()
    else:
        patient6_samples = ['PT6{}_{}'.format(letter, side) for letter in ('A', 'B') for side in ('L', 'R')]
        patient6_mutations = {samplename: all_mut_raw_somatic_filter[all_mut_raw_somatic_filter['sample']==samplename].mutation_id.to_list() for samplename in patient6_samples}
        R_mut = set(patient6_mutations['PT6A_R']).union(set(patient6_mutations['PT6B_R']))
        L_mut = set(patient6_mutations['PT6A_L']).union(set(patient6_mutations['PT6B_L']))
        sample_list = ['PT6_R', 'PT6_L']

    ax = plt.subplot(gs[0, (patient-1)*2:(patient-1)*2+2])
    out = venn2([set(R_mut), set(L_mut)], set_labels=sample_list, ax=ax,
               set_colors=('#961253', '#c54e3a'), alpha=0.5)
    for i, text in enumerate(out.set_labels):
        cur_x, cur_y = text.get_position()
        if i==0:
            text.set_fontsize(15)
            text.set_position([cur_x + 0.2, cur_y])
        else:
            text.set_fontsize(15)
            text.set_position([cur_x - 0.2, cur_y])

    for text in out.subset_labels:
        if text is not None:
            text.set_fontsize(15)


from matplotlib import pyplot, transforms

for ii, side in enumerate(['1_L', '3_L', '3_R', '5_R']):
    sample_list = ['{}{}'.format(sampletype, side) for sampletype in ('PT', 'RD')]
    sub_mut_df = all_mut_raw_somatic_filter[all_mut_raw_somatic_filter['sample'].isin(sample_list)]
    PT_mut = all_mut_raw_somatic_filter[all_mut_raw_somatic_filter['sample']=='PT{}'.format(side)].mutation_id.to_list()
    RD_mut = all_mut_raw_somatic_filter[all_mut_raw_somatic_filter['sample']=='RD{}'.format(side)].mutation_id.to_list()
    patient_nb = int(side[0])
    if 'L' in side:
        patient_side = 1
    else:
        patient_side = 0
    ax = plt.subplot(gs[2, (patient_nb-1)*2+patient_side])
    ax2 = plt.subplot(gs[3, (patient_nb-1)*2+patient_side])
    ax2.set_axis_off()


    out = venn2([set(PT_mut), set(RD_mut)], set_labels=sample_list, ax=ax,
               set_colors=('#a5d060', '#50b1bd'), alpha=0.7)
    t2 = mpl.transforms.Affine2D().rotate_deg(270) + ax.transData
    bottom, top = ax.get_ylim()
    ax.set_ylim(bottom - 0.5, top + 0.5)
    left, right = ax.get_xlim()
    for p in ax.patches:
        p.set_transform(t2)
    for i, text in enumerate(out.set_labels):
        text.set_fontsize(15)
        cur_x, cur_y = text.get_position()
        if i==0:
            text.set_position([(left + right)/2 + 0.3, top+0.3])
        elif i==1:
            text.set_position([(left + right)/2 - 0.3, bottom-0.1])
        
    for i, text in enumerate(out.subset_labels):
        text.set_fontsize(15)
        cur_x, cur_y = text.get_position()
        if i==0:
            text.set_position([(left + right)/2, top-0.15])
        elif i==1:
            text.set_position([(left + right)/2, bottom+0.1])
        elif i==2:
            text.set_position([(left + right)/2, cur_y-0.05])
plt.savefig('figures/fig2_panelB.pdf', bbox_inches="tight")


####################
# figure 2 Panel C #
####################

patient6_samples = ['PT6{}_{}'.format(letter, side) for letter in ('A', 'B') for side in ('L', 'R')]
patient6_mutations = {samplename: all_mut_raw_somatic_filter[all_mut_raw_somatic_filter['sample']==samplename].mutation_id.to_list() for samplename in patient6_samples}


fig = plt.figure(figsize=(10, 6))
gs = gridspec.GridSpec(2, 2)
gs.update(wspace=0.01, hspace=0.0001)

ax = plt.subplot(gs[0, 0:2])
out = venn2([set(patient6_mutations['PT6A_R']).union(set(patient6_mutations['PT6B_R'])),
            set(patient6_mutations['PT6A_L']).union(set(patient6_mutations['PT6B_L']))],
            set_labels=['PT6_R', 'PT6_L'] , ax=ax,
           set_colors=('#961253', '#c54e3a'), alpha=0.5)

for i, text in enumerate(out.set_labels):
    text.set_fontsize(15)

for text in out.subset_labels:
    try:
        text.set_fontsize(15)
    except:
        pass
    
ax = plt.subplot(gs[1, 0])
out = venn2([set(patient6_mutations['PT6A_R']), set(patient6_mutations['PT6B_R'])],
            set_labels=['PT6A_R', 'PT6B_R'] , ax=ax,
           set_colors=('#fdbf60', '#ffde76'), alpha=0.9)
for i, text in enumerate(out.set_labels):
    text.set_fontsize(15)

for text in out.subset_labels:
    try:
        text.set_fontsize(15)
    except:
        pass
    
    
ax = plt.subplot(gs[1, 1])
out = venn2([set(patient6_mutations['PT6A_L']), set(patient6_mutations['PT6B_L'])],
            set_labels=['PT6A_L', 'PT6B_L'] , ax=ax,
           set_colors=('#fdbf60', '#ffde76'), alpha=0.9)
for i, text in enumerate(out.set_labels):
    text.set_fontsize(15)

for text in out.subset_labels:
    try:
        text.set_fontsize(15)
    except:
        pass
    
plt.savefig('figures/fig2_panelC.pdf', bbox_inches="tight")

####################
# figure 3 Panel A #
####################


from ngs_toolbox import *
sns.set_style('white')

chr_p_df = pd.DataFrame.from_dict(centromeres_hg19, orient='index', columns=['centromere'])
chr_p_df = chr_p_df.assign(stop=chr_p_df.centromere*10**6,
                             chromosome=chr_p_df.index.astype(str),
                             start=1,
                            band='p')
chr_q_df = pd.DataFrame.from_dict(centromeres_hg19, orient='index', columns=['centromere'])
chr_q_df = chr_q_df.assign(start=chr_q_df.centromere*10**6,
                            chromosome=chr_q_df.index.astype(str),
                            stop=chromosome_size_hg19.values(),
                           band='q')
cols = ['chromosome', 'start', 'stop', 'band']
chr_pq_df = pd.concat((chr_p_df[cols], chr_q_df[cols])).sort_values(['chromosome', 'start']).reset_index()
chr_pq_df.loc[chr_pq_df.chromosome=="23", 'chromosome'] = 'X'
chr_pq_df.loc[chr_pq_df.chromosome=="24", 'chromosome'] = 'Y'

fig, axes = plt.subplots(nrows=5, figsize=(30, 20))
color_dict = {'gain': '#e84a3f',
              'amplification': '#9c1309',
              'loss': '#4086e3',
              'deletion': '#063a80',
              'CN-LOH': '#d6d6d6'}
vocab_dict = {0: 'deletion', 1: 'loss', 2: 'CN-LOH', 3: 'gain', 4: 'gain', 5: 'amplification'}
for patient_number in range(1, 6):
    cnv_name_1 = 'results/superFreq/plots/P{}_R/data/CNAsegments_PT{}_R.tsv'.format(patient_number, patient_number)
    cnv_name_2 = 'results/superFreq/plots/P{}_L/data/CNAsegments_PT{}_L.tsv'.format(patient_number, patient_number)
    cnv1 = pd.read_csv(cnv_name_1, sep='\t', usecols=['chr', 'start', 'end', 'call', 'clonality', 'clonalityError', 'sigma', 'pCall', 'subclonality', 'subclonalityError'])
    cnv2 = pd.read_csv(cnv_name_2, sep='\t', usecols=['chr', 'start', 'end', 'call', 'clonality', 'clonalityError', 'sigma', 'pCall', 'subclonality', 'subclonalityError'])
    cnv1_band = pd.merge(cnv1, chr_pq_df, left_on='chr', right_on='chromosome', suffixes=['', '_band'])
    cnv1_band = cnv1_band[((cnv1_band.chr==cnv1_band.chromosome)&(cnv1_band.start_band<=cnv1_band.start)&(cnv1_band.start<=cnv1_band.stop))|((cnv1_band.chr==cnv1_band.chromosome)&(cnv1_band.end<=cnv1_band.stop)&(cnv1_band.end>=cnv1_band.start_band))]
    cnv1_band = cnv1_band.assign(start_seg_band=cnv1_band.apply(lambda x: max(x.start, x.start_band), axis=1),
                                 end_seg_band=cnv1_band.apply(lambda x: min(x.end, x.stop), axis=1))

    cnv2_band = pd.merge(cnv2, chr_pq_df, left_on='chr', right_on='chromosome', suffixes=['', '_band'])
    cnv2_band = cnv2_band[((cnv2_band.chr==cnv2_band.chromosome)&(cnv2_band.start_band<=cnv2_band.start)&(cnv2_band.start<=cnv2_band.stop))|((cnv2_band.chr==cnv2_band.chromosome)&(cnv2_band.end<=cnv2_band.stop)&(cnv2_band.end>=cnv2_band.start_band))]
    cnv2_band = cnv2_band.assign(start_seg_band=cnv2_band.apply(lambda x: max(x.start, x.start_band), axis=1),
                                 end_seg_band=cnv2_band.apply(lambda x: min(x.end, x.stop), axis=1))

    cnv1_band = cnv1_band.assign(length=cnv1_band.end_seg_band - cnv1_band.start_seg_band)
    cnv2_band = cnv2_band.assign(length=cnv2_band.end_seg_band - cnv2_band.start_seg_band)
    cnv1_band['new_start'] = cnv1_band.apply(
        lambda row: update_pos_hg19('chr'+row['chr'], row['start_seg_band']),
        axis=1)
    cnv2_band['new_start'] = cnv2_band.apply(
        lambda row: update_pos_hg19('chr'+row['chr'], row['start_seg_band']),
        axis=1)
    cnv1_band = cnv1_band.assign(total_cn=cnv1_band.call.str.replace('?', '').str.replace('CL', '').str.len())
    cnv1_band.loc[cnv1_band.total_cn>5, "total_cn"] = 5
    cnv1_band = cnv1_band[cnv1_band.call!='AB']
    cnv1_band = cnv1_band.assign(name_call=cnv1_band.total_cn.map(vocab_dict))
    cnv1_band = cnv1_band.assign(color=cnv1_band.name_call.map(color_dict))
    cnv1_band_p = cnv1_band[cnv1_band.band=='p']
    cnv1_band_q = cnv1_band[cnv1_band.band=='q']

    cnv2_band = cnv2_band.assign(total_cn=cnv2_band.call.str.replace('?', '').str.replace('CL', '').str.len())
    cnv2_band.loc[cnv2_band.total_cn>5, "total_cn"] = 5
    cnv2_band = cnv2_band[cnv2_band.call!='AB']
    cnv2_band = cnv2_band.assign(name_call=cnv2_band.total_cn.map(vocab_dict))
    cnv2_band = cnv2_band.assign(color=cnv2_band.name_call.map(color_dict))
    cnv2_band_p = cnv2_band[cnv2_band.band=='p']
    cnv2_band_q = cnv2_band[cnv2_band.band=='q']

    ax=axes[5-patient_number]
    
    ax.broken_barh(list(zip(cnv2_band_p.new_start.to_list(), cnv2_band_p.length.to_list())), (20, 9), facecolors=cnv2_band_p.color.to_list())
    ax.broken_barh(list(zip(cnv2_band_q.new_start.to_list(), cnv2_band_q.length.to_list())), (21, 9), facecolors=cnv2_band_q.color.to_list())
    
    ax.broken_barh(list(zip(cnv1_band_p.new_start.to_list(), cnv1_band_p.length.to_list())), (10, 9), facecolors=cnv1_band_p.color.to_list())
    ax.broken_barh(list(zip(cnv1_band_q.new_start.to_list(), cnv1_band_q.length.to_list())), (11, 9), facecolors=cnv1_band_q.color.to_list())

    ax.set_ylim(7, 32)

    ax.set_xlabel('')
    ax.set_ylabel('Patient {}'.format(patient_number), fontsize=35, labelpad=150)
    ax.set_yticks([15, 25])
    ax.tick_params(axis='y', which='both', length=0, pad=50)
    ax.set_yticklabels(['PT{}_R'.format(patient_number), 'PT{}_L'.format(patient_number)], fontsize=25, rotation=0, ha='center')
    plt_genomic_contour(ax, gender='female', lw_factor=2)
    if patient_number==1:
        ax.set_xticklabels([''.format(i) for i in range(1, 25)])
        ax.set_xticklabels(['{}'.format(i) for i in range(1, 23)] + ['X', 'Y'], fontsize=28, rotation=0)
        ax.tick_params(axis='x', which='both', length=0, pad=10)
        
    else:
        ax.set_xticks([])
plt.subplots_adjust(hspace=0.05)
plt.savefig('figures/fig3_panelA.pdf', bbox_inches="tight")


####################
# figure 3 Panel B #
####################


patient6_samples = ['PT6{}_{}'.format(letter, side) for side in ('R', 'L') for letter in ('A', 'B')]
fig, ax = plt.subplots(nrows=1, figsize=(30, 9))
for ii, pn in enumerate(patient6_samples):
    side = pn.split('_')[1]
    cnv_name = 'results/superFreq/plots/P6_{}/data/CNAsegments_{}.tsv'.format(side, pn)
    cnv1 = pd.read_csv(cnv_name, sep='\t', usecols=['chr', 'start', 'end', 'call', 'clonality', 'clonalityError', 'sigma', 'pCall', 'subclonality', 'subclonalityError'])
    cnv1 = cnv1.assign(length=cnv1.end - cnv1.start)
    cnv1['new_start'] = cnv1.apply(
        lambda row: update_pos_hg19('chr'+row['chr'], row['start']),
        axis=1)
    cnv1 = cnv1.assign(total_cn=cnv1.call.str.replace('?', '').str.replace('CL', '').str.len())
    cnv1.loc[cnv1.total_cn>5, "total_cn"] = 5
    cnv1 = cnv1[cnv1.call!='AB']
    cnv1 = cnv1.assign(name_call=cnv1.total_cn.map(vocab_dict))
    cnv1 = cnv1.assign(color=cnv1.name_call.map(color_dict))
    ax.broken_barh(list(zip(cnv1.new_start.to_list(), cnv1.length.to_list())), ((ii+1)*10, 9), facecolors=cnv1.color.to_list())

    ax.set_ylim(7, 52)

ax.set_xlabel('')
ax.set_yticks([15, 25, 35, 45])
ax.set_yticklabels(patient6_samples, fontsize=30)
plt_genomic_contour(ax, gender='female', lw_factor=2)
ax.set_xticklabels([''.format(i) for i in range(1, 25)])
ax.set_xticklabels(['{}'.format(i) for i in range(1, 23)] + ['X', 'Y'], fontsize=28, rotation=0)
ax.tick_params(axis='x', which='both', length=0, pad=10)
plt.savefig('figures/fig3_panelB.pdf', bbox_inches="tight")

####################
# figure 3 Panel C #
####################

fig, axes = plt.subplots(nrows=4, figsize=(30, 16))
sample_list = ['1_L', '3_L', '3_R', '5_R']
for ii, side in enumerate(sample_list):
    cnv_name_1 = 'results/superFreq/plots/P{}/data/CNAsegments_PT{}.tsv'.format(side, side)
    cnv_name_2 = 'results/superFreq/plots/P{}/data/CNAsegments_RD{}.tsv'.format(side, side)
    cnv1 = pd.read_csv(cnv_name_1, sep='\t', usecols=['chr', 'start', 'end', 'call', 'clonality', 'clonalityError', 'sigma', 'pCall', 'subclonality', 'subclonalityError'])
    cnv2 = pd.read_csv(cnv_name_2, sep='\t', usecols=['chr', 'start', 'end', 'call', 'clonality', 'clonalityError', 'sigma', 'pCall', 'subclonality', 'subclonalityError'])
    cnv1 = cnv1.assign(length=cnv1.end - cnv1.start)
    cnv2 = cnv2.assign(length=cnv2.end - cnv2.start)
    cnv1['new_start'] = cnv1.apply(
        lambda row: update_pos_hg19('chr'+row['chr'], row['start']),
        axis=1)
    cnv2['new_start'] = cnv2.apply(
        lambda row: update_pos_hg19('chr'+row['chr'], row['start']),
        axis=1)
    cnv1 = cnv1.assign(total_cn=cnv1.call.str.replace('?', '').str.replace('CL', '').str.len())
    cnv1.loc[cnv1.total_cn>5, "total_cn"] = 5
    cnv1 = cnv1[cnv1.call!='AB']
    cnv1 = cnv1.assign(name_call=cnv1.total_cn.map(vocab_dict))
    cnv1 = cnv1.assign(color=cnv1.name_call.map(color_dict))

    cnv2 = cnv2.assign(total_cn=cnv2.call.str.replace('?', '').str.replace('CL', '').str.len())
    cnv2.loc[cnv2.total_cn>5, "total_cn"] = 5
    cnv2 = cnv2[cnv2.call!='AB']
    cnv2 = cnv2.assign(name_call=cnv2.total_cn.map(vocab_dict))
    cnv2 = cnv2.assign(color=cnv2.name_call.map(color_dict))

    ax=axes[ii]
    ax.broken_barh(list(zip(cnv1.new_start.to_list(), cnv1.length.to_list())), (20, 9), facecolors=cnv1.color.to_list())
    ax.broken_barh(list(zip(cnv2.new_start.to_list(), cnv2.length.to_list())), (10, 9), facecolors=cnv2.color.to_list())

    ax.set_ylim(7, 32)

    ax.set_xlabel('')
    ax.set_yticks([15, 25])
    ax.set_yticklabels(['RD{}'.format(side), 'PT{}'.format(side)], fontsize=30)
    plt_genomic_contour(ax, gender='female', lw_factor=2)
    ax.set_xticklabels([''.format(i) for i in range(1, 25)])
    ax.tick_params(axis='x', which='both', length=0, pad=10)
    if ii < 3:
        ax.set_xticks([])
ax.set_xticklabels(['{}'.format(i) for i in range(1, 23)] + ['X', 'Y'], fontsize=28, rotation=0)
legend_items = list()
for key, color in color_dict.items():
    legend_items.append(Patch(facecolor=color, edgecolor='black', label=key))
ax.legend(handles=legend_items,
          bbox_to_anchor=(0.1, -0.6, 1., .102), fontsize=30, loc=3, ncol=5)
plt.subplots_adjust(hspace=0.05)
plt.savefig('figures/fig3_panelC.pdf', bbox_inches="tight")


#############
# figure S6 #
#############

sns.set_style('white')
all_genes_restr_dedup.index = all_genes_restr_dedup.new_start
cyto = all_genes_restr_dedup.copy()
all_genes_restr_dedup.index = all_genes_restr_dedup.new_end - 1
aa = all_genes_restr_dedup.append(cyto).sort_index()
nb_sample = len(gain_cols)
save_fragl(aa, 'fragl_PT_bilat', nb_sample, colnames=['total_gain', 'total_loss'], format_ext=[])
plt.savefig('figS6.pdf', bbox_inches="tight")


#############
# ext fig 6 #
#############
color_dict = collections.OrderedDict()
color_dict['SBS1']= '#F2786D'
color_dict['SBS2']= '#D86009'
color_dict['SBS13']= '#D86009'
color_dict['SBS3']= '#DA8C05'
color_dict['SBS4']= '#A9A307'
color_dict['SBS5']= '#36B608'
color_dict['SBS6']= '#07BC7C'
color_dict['SBS7a']= '#f2b67c'
color_dict['SBS7b']= '#f2b67c'
color_dict['SBS7c']= '#f2b67c'
color_dict['SBS7d']= '#f2b67c'
color_dict['SBS8']= '#02BFC0'
color_dict['SBS9']= '#00AEF7'
color_dict['SBS10a']= '#66A428'
color_dict['SBS10b']= '#66A428'
color_dict['SBS11']= '#938CFF'
color_dict['SBS12']= '#E071EC'
color_dict['SBS14']= '#F566BE'
color_dict['SBS15']= '#7CCA7C'
color_dict['SBS16']= '#356AB2'
color_dict['SBS17a']= '#C15B00'
color_dict['SBS17b']= '#C15B00'
color_dict['SBS18']= '#666666'
color_dict['SBS19']= '#756EB5'
color_dict['SBS20']= '#E7AC00'
color_dict['SBS21']= '#A5CEE4'
color_dict['SBS22']= '#2CA121'
color_dict['SBS23']= '#335a9c'
color_dict['SBS24']= '#FEC068'
color_dict['SBS25']= '#a5f0b1'
color_dict['SBS26']= '#6A399C'
color_dict['SBS28']= '#E71408'
color_dict['SBS29']= '#4e88a6'
color_dict['SBS30']= '#994BA5'
color_dict['SBS31']= '#0cf531'
color_dict['SBS33']= '#A75620'
color_dict['SBS34']= '#62C3A4'
color_dict['SBS35']= '#E988C4'
color_dict['SBS36']= '#E6C591'
color_dict['SBS37']= '#BEADD5'
color_dict['SBS38']= '#F3007E'
color_dict['SBS39']= '#089F76'
color_dict['SBS40']= '#A77709'
color_dict['SBS41']= '#a80526'
color_dict['SBS42']= 'black'
color_dict['SBS43']= '#6b088a'
color_dict['SBS44']= '#fce428'
color_dict['SBS45']= '#f0d0f2'
color_dict['SBS46']= '#78aaff'
color_dict['SBS47']= '#ffb108'
color_dict['SBS49']= '#c108ff'
color_dict['SBS32']= '#c27289'
color_dict['SBS54']= '#1320d1'
color_dict['SBS56']= '#e8dba2'
color_dict['SBS58']= '#b0eb31'
color_dict['SBS59']= '#b8b8b8'

used_signatures = list()

sig_table = pd.read_csv('bilat_decompTumor2Sig.csv', sep='\t')

fig = plt.figure(constrained_layout=False, figsize=(25, 18))
gs = fig.add_gridspec(ncols=6, nrows=2, wspace=0.2)


for patient_number in range(1, 6):
    sig_df = sig_table[['PT{}_R'.format(patient_number), 'PT{}_L'.format(patient_number)]].T
    sig_cols = ListedColormap([color_dict[s] for s in sig_df.columns.to_list()])
    ax = fig.add_subplot(gs[0, patient_number-1])
    sig_df.plot.bar(stacked=True, cmap=sig_cols, ax=ax, rot=0, width=0.7)
    ax.set_ylim([0, 1])
    ax.get_legend().remove()
    if patient_number != 1:
        ax.get_yaxis().set_visible(False)
    else:
        ax.set_ylabel('signature proportion', fontsize=30)

legend_items = list()
for key in ('SBS1', 'SBS2', 'SBS3', 'SBS5', 'SBS8', 'SBS9', 'SBS13', 'SBS17a', 'SBS17b', 'SBS18', 'SBS37', 'SBS40', 'SBS41'):
    legend_items.append(Patch(facecolor=color_dict[key], label=key))



for ii, side in enumerate(['1_L', '3_L', '3_R', '5_R']):
    sig_df = sig_table[['PT{}'.format(side), 'RD{}'.format(side)]].T
    sig_cols = ListedColormap([color_dict[s] for s in sig_df.columns.to_list()])
    ax = fig.add_subplot(gs[1, ii])
    sig_df.plot.bar(stacked=True, cmap=sig_cols, ax=ax, rot=0, width=0.7)
    ax.set_ylim([0, 1])
    used_signatures += sig_df.columns[sig_df.T.sum(axis=1)>0].to_list()
    ax.get_legend().remove()


    if ii!=0:
        ax.get_yaxis().set_visible(False)
    else:
        ax.set_ylabel('signature proportion', fontsize=30)
        ax.legend(handles=legend_items,
          bbox_to_anchor=(-0.6, -0.35, 1., .102),
          fontsize=25, loc=3, ncol=7)
        

plt.savefig('figures/figS6.pdf', bbox_inches="tight")

used_signatures = list()
patient6_samples = ['PT6{}_{}'.format(letter, side) for side in ('R', 'L') for letter in ('A', 'B')]
fig, ax = plt.subplots(figsize=(8, 7))

sig_df = sig_table[patient6_samples].T

sig_cols = ListedColormap([color_dict[s] for s in sig_df.columns.to_list()])
sig_df.plot.bar(stacked=True, cmap=sig_cols, ax=ax, rot=0, width=0.6)
ax.set_ylim([0, 1])
used_signatures += sig_df.columns[sig_df.T.sum(axis=1)>0].to_list()
ax.get_legend().remove()
ax.set_ylabel('signature proportion', fontsize=30)

legend_items = list()
for key in ('SBS1', 'SBS2', 'SBS3', 'SBS5', 'SBS8', 'SBS9', 'SBS13', 'SBS17a', 'SBS17b', 'SBS18', 'SBS37', 'SBS40', 'SBS41'):
    legend_items.append(Patch(facecolor=color_dict[key], label=key))

ax.legend(handles=legend_items,
  bbox_to_anchor=(-0.2, -0.7, 1., .102),
  fontsize=25, loc=3, ncol=3)
plt.savefig('figures/figS6_panelC.pdf', bbox_inches="tight")



