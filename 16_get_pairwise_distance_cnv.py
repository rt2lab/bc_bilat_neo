from sklearn.metrics.pairwise import *

sample_list = ["PT1_R", "PT1_L", "RD1_L",
             "PT2_R", "PT2_L",
             "PT3_R", "RD3_R", "PT3_L", "RD3_L",
             "PT4_R", "PT4_L",
             "PT5_R", "RD5_R", "PT5_L",
             "PT6A_R", "PT6B_R", "PT6A_L", "PT6B_L"]

cnv_dict = dict()

for sample_name in sample_list:
    cnv_name = 'data/pre_fragl_{}.tsv'.format(patientside, sample_name)
    cnv_table = pd.read_csv(cnv_name, sep='\t')
    cnv_dict[sample_name] = cnv_table.avg_tot_cn.values

big_cnv_df = pd.DataFrame.from_dict(cnv_dict, orient='index')

cos_df = pd.DataFrame(cosine_similarity(big_cnv_df.values-2),
                      index=sample_list, columns=sample_list)

sns.set_context('poster', font_scale=0.6)
fig, ax = plt.subplots(figsize=(17,15))
sns.heatmap(cos_df, annot=True, fmt='.2f', ax=ax, center=0, cmap="vlag", vmin=-1, vmax=1)
plt.savefig('pairwise_cosine_similarity_original_cohort.pdf', bbox_inches="tight")