#!/usr/bin/python
import subprocess, sys, getopt
from collections import defaultdict
from collections import OrderedDict
from statistics import mean, median
from itertools import groupby
from operator import itemgetter
from collections import Counter
import re


def classify(val):
    res = []
    if val > 0:
        res.append('n_bp_0X')
    if val > 1:
        res.append('n_bp_1X')
    if val > 2:
        res.append('n_bp_2X')
    if val > 3:
        res.append('n_bp_3X')
    if val > 4:
        res.append('n_bp_4X')
    if val > 5:
        res.append('n_bp_5X')
    if val > 10:
        res.append('n_bp_10X')
    if val > 20:
        res.append('n_bp_20X')
    if val > 30:
        res.append('n_bp_30X')
    if val > 40:
        res.append('n_bp_40X')
    if val > 50:
        res.append('n_bp_50X')
    if val > 60:
        res.append('n_bp_60X')
    if val > 70:
        res.append('n_bp_70X')
    if val > 80:
        res.append('n_bp_80X')
    if val > 90:
        res.append('n_bp_90X')
    if val > 100:
        res.append('n_bp_100X')
    return res

def generateDepth(out_file_path):

    depthByBase = defaultdict(lambda: defaultdict(int))
    #cmd = f"samtools depth -a -b {bed_file_path} -q 20 -Q 20 -s {bam_file_path}"
    #p = subprocess.Popen(cmd, shell= True,stdout=subprocess.PIPE)
    #for line in iter(p.stdout.readline, b''):
    summaryRow = defaultdict(int)
    for line in sys.stdin:
        words = line.split()
        print(words)
        if len(words) == 3:
            depthByBase[words[0]][int(words[1])] = int(words[2]) 

    with open(out_file_path,"w") as out:
        out.write("CHR"+"\t"+"region_start"+"\t"+"region_end"+"\t"+"region_length"+"\t"+"n_bp_0X"+"\t"+"n_bp_1X"+"\t"+"n_bp_2X"+"\t"+"n_bp_3X"+"\t"+"n_bp_4X"+"\t"
                    +"n_bp_5X"+"\t"+"n_bp_10X"+"\t"+"n_bp_20X"+"\t"+"n_bp_30X"+"\t"+"n_bp_40X"+"\t"+"n_bp_50X"+"\t"+"n_bp_60X"+"\t"+"n_bp_70X"+"\t"+"n_bp_80X"+"\t"
                    +"n_bp_90X"+"\t"+"n_bp_100X"+"\t"+"average_depth"+"\t"+"sum_depth"+"\n")
        totalBases = 0
        for key in depthByBase.keys():
            for k, g in groupby(enumerate(depthByBase[key]), key=lambda x:x[0]-x[1]):
                coords = list(map(itemgetter(1),g))
                res = dict((k, depthByBase[key][k]) for k in coords)
                classes = Counter(classification for val in res.values() for classification in classify(val))
                totalBases += coords[-1]-coords[0]+1
                summaryRow['n_bp_0X'] += classes['n_bp_0X']
                summaryRow['n_bp_1X'] += classes['n_bp_1X']
                summaryRow['n_bp_2X'] += classes['n_bp_2X']
                summaryRow['n_bp_3X'] += classes['n_bp_3X']
                summaryRow['n_bp_4X'] += classes['n_bp_4X']
                summaryRow['n_bp_5X'] += classes['n_bp_5X']
                summaryRow['n_bp_10X'] += classes['n_bp_10X']
                summaryRow['n_bp_20X'] += classes['n_bp_20X']
                summaryRow['n_bp_30X'] += classes['n_bp_30X']
                summaryRow['n_bp_40X'] += classes['n_bp_40X']
                summaryRow['n_bp_50X'] += classes['n_bp_50X']
                summaryRow['n_bp_60X'] += classes['n_bp_60X']
                summaryRow['n_bp_70X'] += classes['n_bp_70X']
                summaryRow['n_bp_80X'] += classes['n_bp_80X']
                summaryRow['n_bp_90X'] += classes['n_bp_90X']
                summaryRow['n_bp_100X'] += classes['n_bp_100X']
                summaryRow['average_depth'] += round(mean(res.values()),2)
                summaryRow['sum_depth'] += sum(res.values())
                out.write(key+"\t"+str(coords[0])+"\t"+str(coords[-1])+"\t"+str(coords[-1]-coords[0]+1)+"\t"+str(classes['n_bp_0X'])+"\t"
                    +str(classes['n_bp_1X'])+"\t"+str(classes['n_bp_2X'])+"\t"+str(classes['n_bp_3X'])+"\t"+str(classes['n_bp_4X'])+"\t"
                    +str(classes['n_bp_5X'])+"\t"+str(classes['n_bp_10X'])+"\t"+str(classes['n_bp_20X'])+"\t"+str(classes['n_bp_30X'])+"\t"+str(classes['n_bp_40X'])+"\t"
                    +str(classes['n_bp_50X'])+"\t"+str(classes['n_bp_60X'])+"\t"+str(classes['n_bp_70X'])+"\t"+str(classes['n_bp_80X'])+"\t"
                    +str(classes['n_bp_90X'])+"\t"+str(classes['n_bp_100X'])+"\t"+str(round(mean(res.values()),2))+"\t"+str(sum(res.values()))+"\n")
        out.write("AUTOSOMAL"+"\t"+"NA"+"\t"+"NA"+"\t"+str(totalBases)+"\t"+str(summaryRow['n_bp_0X'])+"\t"
                +str(summaryRow['n_bp_1X'])+"\t"+str(summaryRow['n_bp_2X'])+"\t"+str(summaryRow['n_bp_3X'])+"\t"+str(summaryRow['n_bp_4X'])+"\t"
                +str(summaryRow['n_bp_5X'])+"\t"+str(summaryRow['n_bp_10X'])+"\t"+str(summaryRow['n_bp_20X'])+"\t"+str(summaryRow['n_bp_30X'])+"\t"+str(summaryRow['n_bp_40X'])+"\t"
                +str(summaryRow['n_bp_50X'])+"\t"+str(summaryRow['n_bp_60X'])+"\t"+str(summaryRow['n_bp_70X'])+"\t"+str(summaryRow['n_bp_80X'])+"\t"
                +str(summaryRow['n_bp_90X'])+"\t"+str(summaryRow['n_bp_100X'])+"\t"+str(round((summaryRow['sum_depth']/totalBases),2))+"\t"+str(summaryRow['sum_depth'])+"\n")


def generateAllDepth(out_file_path):

    outStats = defaultdict(lambda: defaultdict(int))
    #depthsForChr = defaultdict(lambda: defaultdict(int))

    #cmd = f"samtools depth -a -q 20 -Q 20 -s {bam_file_path}"
    #p = subprocess.Popen(cmd, shell= True,stdout=subprocess.PIPE)
    #for line in iter(p.stdout.readline, b''):
    summaryRow = defaultdict(int)
    totalBases = 0
    for line in sys.stdin:
            words = line.split()
            if len(words) == 3:
                if outStats[words[0]]['region_start'] == 0:
                    outStats[words[0]]['region_start'] = int(words[1])
                outStats[words[0]]['region_length']+=1
                depth = int(words[2])
                #depthsForChr[words[0]][depth]+=1
                outStats[words[0]]['sum_depth']+=depth
                if depth > 0: 
                    outStats[words[0]]['n_bp_0X']+= 1
                if depth > 1: 
                    outStats[words[0]]['n_bp_1X']+= 1
                if depth > 2: 
                    outStats[words[0]]['n_bp_2X']+= 1
                if depth > 3: 
                    outStats[words[0]]['n_bp_3X']+= 1
                if depth > 4: 
                    outStats[words[0]]['n_bp_4X']+= 1
                if depth > 5: 
                    outStats[words[0]]['n_bp_5X']+= 1
                if depth > 10: 
                    outStats[words[0]]['n_bp_10X']+= 1
                if depth > 20: 
                    outStats[words[0]]['n_bp_20X']+= 1
                if depth > 30: 
                    outStats[words[0]]['n_bp_30X']+= 1
                if depth > 40: 
                    outStats[words[0]]['n_bp_40X']+= 1
                if depth > 50: 
                    outStats[words[0]]['n_bp_50X']+= 1
                if depth > 60: 
                    outStats[words[0]]['n_bp_60X']+= 1
                if depth > 70: 
                    outStats[words[0]]['n_bp_70X']+= 1
                if depth > 80: 
                    outStats[words[0]]['n_bp_80X']+= 1
                if depth > 90: 
                    outStats[words[0]]['n_bp_90X']+= 1
                if depth > 100: 
                    outStats[words[0]]['n_bp_100X']+= 1
            else:
                outStats[words[0]]['region_length']+=1

    with open(out_file_path,"w") as out:
        out.write("CHR"+"\t"+"region_start"+"\t"+"region_end"+"\t"+"region_length"+"\t"+"n_bp_0X"+"\t"+"n_bp_1X"+"\t"+"n_bp_2X"+"\t"+"n_bp_3X"+"\t"+"n_bp_4X"+"\t"
                    +"n_bp_5X"+"\t"+"n_bp_10X"+"\t"+"n_bp_20X"+"\t"+"n_bp_30X"+"\t"+"n_bp_40X"+"\t"+"n_bp_50X"+"\t"+"n_bp_60X"+"\t"+"n_bp_70X"+"\t"+"n_bp_80X"+"\t"
                    +"n_bp_90X"+"\t"+"n_bp_100X"+"\t"+"average_depth"+"\t"+"sum_depth"+"\n")

        for key in outStats.keys():
            if (re.match("chr[1-9][0-2]?$",key)):
                summaryRow['length'] += outStats[key]["region_length"]
                summaryRow['n_bp_0X'] += outStats[key]['n_bp_0X']
                summaryRow['n_bp_1X'] += outStats[key]['n_bp_1X']
                summaryRow['n_bp_2X'] += outStats[key]['n_bp_2X']
                summaryRow['n_bp_3X'] += outStats[key]['n_bp_3X']
                summaryRow['n_bp_4X'] += outStats[key]['n_bp_4X']
                summaryRow['n_bp_5X'] += outStats[key]['n_bp_5X']
                summaryRow['n_bp_10X'] += outStats[key]['n_bp_10X']
                summaryRow['n_bp_20X'] += outStats[key]['n_bp_20X']
                summaryRow['n_bp_30X'] += outStats[key]['n_bp_30X']
                summaryRow['n_bp_40X'] += outStats[key]['n_bp_40X']
                summaryRow['n_bp_50X'] += outStats[key]['n_bp_50X']
                summaryRow['n_bp_60X'] += outStats[key]['n_bp_60X']
                summaryRow['n_bp_70X'] += outStats[key]['n_bp_70X']
                summaryRow['n_bp_80X'] += outStats[key]['n_bp_80X']
                summaryRow['n_bp_90X'] += outStats[key]['n_bp_90X']
                summaryRow['n_bp_100X'] += outStats[key]['n_bp_100X']
                outStats[key]['region_end'] = outStats[key]['region_start']+outStats[key]['region_length']-1
                outStats[key]['average_depth'] = round((outStats[words[0]]['sum_depth']/outStats[key]['region_length']),2)
                summaryRow['average_depth'] += outStats[key]['average_depth']
                summaryRow['sum_depth'] += outStats[key]['sum_depth']
                out.write(key+"\t"+str(outStats[key]["region_start"])+"\t"+str(outStats[key]["region_end"])+"\t"+str(outStats[key]["region_length"])+"\t"+str(outStats[key]["n_bp_0X"])+"\t"
                    +str(outStats[key]["n_bp_1X"])+"\t"+str(outStats[key]["n_bp_2X"])+"\t"+str(outStats[key]["n_bp_3X"])+"\t"+str(outStats[key]["n_bp_4X"])+"\t"
                    +str(outStats[key]["n_bp_5X"])+"\t"+str(outStats[key]["n_bp_10X"])+"\t"+str(outStats[key]["n_bp_20X"])+"\t"+str(outStats[key]["n_bp_30X"])+"\t"+str(outStats[key]["n_bp_40X"])+"\t"
                    +str(outStats[key]["n_bp_50X"])+"\t"+str(outStats[key]["n_bp_60X"])+"\t"+str(outStats[key]["n_bp_70X"])+"\t"+str(outStats[key]["n_bp_80X"])+"\t"
                    +str(outStats[key]["n_bp_90X"])+"\t"+str(outStats[key]["n_bp_100X"])+"\t"+str(outStats[key]["average_depth"])+"\t"+str(outStats[words[0]]['sum_depth'])+"\n")
        out.write("AUTOSOMAL"+"\t"+"NA"+"\t"+"NA"+"\t"+str(summaryRow['length'])+"\t"+str(summaryRow['n_bp_0X'])+"\t"
                +str(summaryRow['n_bp_1X'])+"\t"+str(summaryRow['n_bp_2X'])+"\t"+str(summaryRow['n_bp_3X'])+"\t"+str(summaryRow['n_bp_4X'])+"\t"
                +str(summaryRow['n_bp_5X'])+"\t"+str(summaryRow['n_bp_10X'])+"\t"+str(summaryRow['n_bp_20X'])+"\t"+str(summaryRow['n_bp_30X'])+"\t"+str(summaryRow['n_bp_40X'])+"\t"
                +str(summaryRow['n_bp_50X'])+"\t"+str(summaryRow['n_bp_60X'])+"\t"+str(summaryRow['n_bp_70X'])+"\t"+str(summaryRow['n_bp_80X'])+"\t"
                +str(summaryRow['n_bp_90X'])+"\t"+str(summaryRow['n_bp_100X'])+"\t"+str(round((summaryRow['sum_depth']/summaryRow['length']),2))+"\t"+str(summaryRow['sum_depth'])+"\n")




def main(argv):
    bam_file_path = ''
    bed_file_path = ''
    out_file_path = ''
    all_autosomal = False
    try:
        opts, args = getopt.getopt(argv,"ao:")
    except getopt.GetoptError:
        sys.exit(2)
    if not opts:
        sys.exit(2)
    for opt, arg in opts:
        if (opt == "-o"):
            out_file_path = arg
        elif (opt == "-a"):
            all_autosomal = True
    if (all_autosomal):
        generateAllDepth(out_file_path)
    else:
        generateDepth(out_file_path)


if __name__ == '__main__':
    main(sys.argv[1:])