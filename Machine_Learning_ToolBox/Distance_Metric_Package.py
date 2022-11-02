import sklearn
import numpy as np
import scipy
from scipy.spatial import distance
from scipy.special import rel_entr
import rdkit

import pandas as pd

# Distance
# 1.Euclidean Distance
def ed(S1, S2):
    return distance.euclidean(S1, S2)
# 2. Manhattan Distance
def md(S1, S2):
    return distance.cityblock(S1, S2)
# 3. Minkowski Distance
def mind(S1, S2):
    return distance.minkowski(S1, S2)
# 4. Hamming Distance
def hd(S1, S2):
    return distance.hamming(S1, S2)
# 5. L1, L2 and Lp norms
# 6. Cosine Distance and Similarity
def cosd(S1, S2):
    return distance.cosine(S1, S2)
# 7. Jaccard
def jdis(S1, S2):
    return distance.jaccard(S1, S2)
# 8. Bray-Curtis
def bcdis(S1, S2):
    return distance.braycurtis(S1, S2)
# 9. Canberra
def cand(S1, S2):
    return distance.canberra(S1, S2)
# 10 Chebyshev
def ched(S1, S2):
    return distance.chebyshev(S1, S2)
# 11 JensenShannon
def jsd(P1, P2):
    return distance.jensenshannon(P1, P2)
#12 Mahalanobis
# def mahd(S1, S2):
#     return distance.mahalanobis(S1, S2)
#13 Seuclidean
# def sed(S1, S2):
#     return distance.seuclidean(S1, S2)
#14 Sqeuclidean
def sqd(S1, S2):
    return distance.sqeuclidean(S1, S2)
#15 Dice
def ddis(S1, S2):
    return distance.dice(S1, S2)
#16 Kulsinski
def kuld(S1, S2):
    return distance.kulsinski(S1, S2)
#17 Rogers-Tanimoto
def rtd(S1, S2):
    return distance.rogerstanimoto(S1, S2)
#18 Russell-Rao
def rrd(S1, S2):
    return distance.russellrao(S1, S2)
#19 Sokal-Michener
def smd(S1, S2):
    return distance.sokalmichener(S1, S2)
#20 Sokal-Sneath
def ssd(S1, S2):
    return distance.sokalsneath(S1, S2)
#21 Yule
def yd(S1, S2):
    return distance.yule(S1, S2)



# 3. Kullback–Leibler
def kldiv(P, Q):
    return sum(rel_entr(P,Q))
# 4. Unifrac
def udis(S1, S2):
    pass
# 5. Bregman Distance
# 6. itacoro-sito distance
# 7. Tanimoto Similarity
def tsim(V1, V2):
    return rdkit.Chem.DataStructs.cDataStructs.TanimotoSimilarity(V1, V2)
# 8. Log-Likelihood Ratios
# 9. Spearman Correlation
# 10. Pearson Correlation
