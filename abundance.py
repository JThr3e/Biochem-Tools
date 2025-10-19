import matplotlib.pyplot as plt
from collections import defaultdict

phosphoric = {
    "name": "phosphoric",
    "species": ["H3PO4", "H2PO4_1-", "HPO4_2-", "PO4_3-"],
    "pKa" : [2.14,7.2, 12.37] 
}

histidine = {
    "name": "Histidine",
    "species": ["Carboxyl", "Carboxylate", "Imidazole", "Amino"],
    "pKa" : [1.8, 6.0, 9.17] 
}



#given a list of pka and pH, return a list of ratio
def henderson_h(pka_list, pH):
    ratio = []
    for pKa in pka_list:
        ratio.append(10**(pH - pKa))
    return ratio

# Normalize the list of ratios to the equivalent fractions
# Essentially does the [A-]/[HA] => [A-]/([HA]+[A-]) & [HA]/([HA]+[A-]) step
def normalize_to_frac(ratio_list):
    frac_list = []
    for i in range(1, len(ratio_list)):
        ratio_list[i] = ratio_list[i]*ratio_list[i-1]
    denom = sum(ratio_list) + 1
    for ratio in ratio_list:
        frac_list.append(ratio/denom)
    frac_list.insert(0, 1.0/denom)
    return frac_list

def gen_chart(acid_desc):
    import numpy as np
    pH = np.linspace(0, 14, 200)
    answer_dict = defaultdict(list)
    for elem in pH:
        abundance = normalize_to_frac(henderson_h(acid_desc["pKa"], elem))
        assert len(abundance) == len(acid_desc["species"])
        for sp,ab in zip(acid_desc["species"], abundance):
            answer_dict[sp].append(ab)
    for sp in acid_desc["species"]:
        plt.plot(pH, answer_dict[sp], label=sp)
    plt.title(acid_desc["name"])
    plt.legend()
    plt.savefig(acid_desc["name"]+'.png') 
    
gen_chart(histidine)            

