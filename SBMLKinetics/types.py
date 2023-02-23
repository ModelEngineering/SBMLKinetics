# This script was written by Jin Xu and available on Github
# https://github.com/ModelEngineering/SBMLKinetics
# This file defines the objects of R_type and K_type.


class R_type:
    def __init__(self, rct_num, prd_num):
        if type(rct_num) == int and type(prd_num) == int and 0<=rct_num<=3 and 0<=prd_num<=3:
            self.rct_num = rct_num
            self.prd_num = prd_num
        else:
            ValueError('Please enter a valid integer number of reactant or product!')

# print(type(R_type(0, 0).rct_num))

class K_type:
    def __init__(self, K_type_str):
        if type(K_type_str) == str and K_type_str in ["ZERO", "UNDR", "UNMO", "BIDR", "BIMO", "MM", \
            "MMCAT", "HILL", "FR", "NA"]:
            self.K_type_str = K_type_str
        else:
            ValueError('Please enter a valid kinetic law type string!')

# print(type(K_type("ZERO").K_type_str))