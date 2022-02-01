from sklearn.metrics import adjusted_mutual_info_score,adjusted_rand_score
from scipy.stats import entropy
from statistics import mean
from statistics import stdev
import random
import scipy.stats as stats
import pandas as pd

def calc_MI(category_1,category_2):
    return adjusted_mutual_info_score(category_1,category_2,average_method='arithmetic')

def calc_fisherExact(contingency_table):
    return stats.fisher_exact(contingency_table)

def calc_shanon_entropy(value_list):
    return

def calc_adjusted_rand(category_1,category_2):
    return adjusted_rand_score(category_1,category_2)


df = pd.read_csv('/Users/jrobertson/Google Drive/__MOB-CLuster_Manuscript/ClusterComparison_MOBsuite_db_versions.txt',sep="\t",header=0)
nmi = calc_adjusted_rand(df['primary_cluster_id'],
              df['mob_v2_cluster'])
print(nmi)
nmi = calc_adjusted_rand(df['primary_cluster_id'],
              df['mob_v1_cluster_single'])
print(nmi)
nmi = calc_adjusted_rand(df['primary_cluster_id'],
              df['mob_v1_cluster_complete'])
print(nmi)