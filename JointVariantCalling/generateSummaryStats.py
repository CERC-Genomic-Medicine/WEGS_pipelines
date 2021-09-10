#!/usr/bin/python
import subprocess, sys, getopt
from collections import defaultdict
from collections import OrderedDict
from statistics import mean, median
from itertools import groupby
from operator import itemgetter
from collections import Counter
import re

def generateStats():

    known = defaultdict(int)
    novel = defaultdict(int)
    for line in sys.stdin:
        words = line.split()
        info_fields = words[7].split(';')
        AC = info_fields[0].split('=')

        if words[2] != '.':
            if words[6] == 'PASS':
                known["PASS"] += 1
                if int(AC[1]) == 1:
                    known["AC_eq_1"] += 1
                elif int(AC[1]) == 2:
                    known["AC_eq_2"] += 1
                elif int(AC[1]) >= 3 and int(AC[1]) <= 5:
                    known["AC_3to5"] += 1
                elif int(AC[1]) >= 6 and int(AC[1]) <= 10:
                    known["AC_6to10"] += 1 
                elif int(AC[1]) >= 11 and int(AC[1]) <= 100:
                    known["AC_11to100"] += 1
                elif int(AC[1]) > 100:
                    known["AC_gt_100"] += 1 
            else:
                known["FAIL"] += 1          
        else:
            if words[6] == 'PASS':
                novel["PASS"] += 1
                if int(AC[1]) == 1:
                    novel["AC_eq_1"] += 1
                elif int(AC[1]) == 2:
                    novel["AC_eq_2"] += 1
                elif int(AC[1]) >= 3 and int(AC[1]) <= 5:
                    novel["AC_3to5"] += 1
                elif int(AC[1]) >= 6 and int(AC[1]) <= 10:
                    novel["AC_6to10"] += 1 
                elif int(AC[1]) >= 11 and int(AC[1]) <= 100:
                    novel["AC_11to100"] += 1
                elif int(AC[1]) > 100:
                    novel["AC_gt_100"] += 1
            else:
                novel["FAIL"] += 1
    
    var_total = known["PASS"]+known["FAIL"]+novel["PASS"]+novel["FAIL"]
    var_pass = known["PASS"]+novel["PASS"]
    var_fail = known["FAIL"]+novel["FAIL"]
    var_known_tot = known["PASS"]+known["FAIL"]
    var_novel_tot = novel["PASS"]+novel["FAIL"]

    print(var_total,',',var_pass,',',var_fail)
    print(var_known_tot,',',known["PASS"],',',known["FAIL"])
    print(var_novel_tot,',',novel["PASS"],',',novel["FAIL"])

    print("\n\n")

    print(known["AC_eq_1"]+novel["AC_eq_1"],',',novel["AC_eq_1"])
    print(known["AC_eq_2"]+novel["AC_eq_2"],',',novel["AC_eq_2"])
    print(known["AC_3to5"]+novel["AC_3to5"],',',novel["AC_3to5"])
    print(known["AC_6to10"]+novel["AC_6to10"],',',novel["AC_6to10"])
    print(known["AC_11to100"]+novel["AC_11to100"],',',novel["AC_11to100"])
    print(known["AC_gt_100"]+novel["AC_gt_100"],',',novel["AC_gt_100"])




    #print('Known Variants')
    #for key, val in known.items():
    #    print(key,': ',val)
    #print('Novel Variants')
    #for key, val in novel.items():
    #    print(key,': ',val)

if __name__ == '__main__':
    generateStats()
        
    