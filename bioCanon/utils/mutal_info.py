from sklearn.metrics import adjusted_mutual_info_score
from scipy.stats import entropy
from statistics import mean
from statistics import stdev
import random
import scipy.stats as stats

def calc_MI(category_1,category_2):
    return adjusted_mutual_info_score(category_1,category_2,average_method='arithmetic')

def calc_fisherExact(contingency_table):
    return stats.fisher_exact(contingency_table)

def calc_shanon_entropy(value_list):
    return

def calc_adjusted_rand():
    return


nmi = calc_MI((29,29,29,29,29,29,29,29,29,29,29,29,29,29,29,29,29,29,29,29,29,29,29),
              (32,32,32,32,32,32,32,32,32,32,32,32,32,32,32,32,32,32,32,32,29,29,29))
print(nmi)