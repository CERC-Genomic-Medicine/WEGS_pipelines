#!/usr/bin/python
import subprocess, sys, getopt
from collections import defaultdict
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
    elif val > 40:
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

def generateDepth(bam_file_path,bed_file_path,out_file_path):

    depthByBase = defaultdict(lambda: defaultdict(int))
    cmd = f"samtools depth -a -b {bed_file_path} -q 20 -Q 20 -s {bam_file_path}"
    p = subprocess.Popen(cmd, shell= True,stdout=subprocess.PIPE)
    i = 0
    l = 0
    for line in iter(p.stdout.readline, b''):
        words = line.decode().split()
        depthByBase[words[0]][int(words[1])] = int(words[2])

    with open(out_file_path,"w") as out:
        out.write("CHR"+"\t"+"region_start"+"\t"+"region_end"+"\t"+"region_length"+"\t"+"n_bp_0X"+"\t"+"n_bp_1X"+"\t"+"n_bp_2X"+"\t"+"n_bp_3X"+"\t"+"n_bp_4X"+"\t"
                    +"n_bp_5X"+"\t"+"n_bp_10X"+"\t"+"n_bp_20X"+"\t"+"n_bp_30X"+"\t"+"n_bp_40X"+"\t"+"n_bp_50X"+"\t"+"n_bp_60X"+"\t"+"n_bp_70X"+"\t"+"n_bp_80X"+"\t"
                    +"n_bp_90X"+"\t"+"n_bp_100X"+"\t"+"average_depth"+"\t"+"median_depth"+"\n")
   
        for key in depthByBase.keys():
            for k, g in groupby(enumerate(depthByBase[key]), key=lambda x:x[0]-x[1]):
                coords = list(map(itemgetter(1),g))
                res = dict((k, depthByBase[key][k]) for k in coords)
                classes = Counter(classification for val in res.values() for classification in classify(val))
                out.write(key+"\t"+str(coords[0])+"\t"+str(coords[-1])+"\t"+str(coords[-1]-coords[0]+1)+"\t"+str(classes['n_bp_0X'])+"\t"
                    +str(classes['n_bp_1X'])+"\t"+str(classes['n_bp_2X'])+"\t"+str(classes['n_bp_3X'])+"\t"+str(classes['n_bp_4X'])+"\t"
                    +str(classes['n_bp_5X'])+"\t"+str(classes['n_bp_10X'])+"\t"+str(classes['n_bp_20X'])+"\t"+str(classes['n_bp_30X'])+"\t"+str(classes['n_bp_40X'])+"\t"
                    +str(classes['n_bp_50X'])+"\t"+str(classes['n_bp_60X'])+"\t"+str(classes['n_bp_70X'])+"\t"+str(classes['n_bp_80X'])+"\t"
                    +str(classes['n_bp_90X'])+"\t"+str(classes['n_bp_100X'])+"\t"+str(round(mean(res.values()),2))+"\t"+str(median(res.values()))+"\n")



def generateAllDepth(bam_file_path,out_file_path):

    outStats = defaultdict(dict)
    coords = defaultdict(list)
    depthsForChr = defaultdict(list)

    cmd = f"samtools depth -a -q 20 -Q 20 -s {bam_file_path}"
    p = subprocess.Popen(cmd, shell= True,stdout=subprocess.PIPE)

    for chr in range(1,23):
        fields = {}
        fields['region_start'] = 0
        fields['region_end'] = 0
        fields['region_length'] = 0
        fields['n_bp_0X'] = 0
        fields['n_bp_1X'] = 0
        fields['n_bp_2X'] = 0
        fields['n_bp_3X'] = 0
        fields['n_bp_4X'] = 0
        fields['n_bp_5X'] = 0
        fields['n_bp_10X'] = 0
        fields['n_bp_20X'] = 0
        fields['n_bp_30X'] = 0
        fields['n_bp_40X'] = 0
        fields['n_bp_50X'] = 0
        fields['n_bp_60X'] = 0
        fields['n_bp_70X'] = 0
        fields['n_bp_80X'] = 0
        fields['n_bp_90X'] = 0
        fields['n_bp_100X'] = 0
        fields['average_depth'] = 0.0
        fields['median_depth'] = 0.0
        outStats['chr'+str(chr)] = fields

    for line in iter(p.stdout.readline, b''):
            words = line.decode().split()
            coords[words[0]].append(int(words[1]))
            depthsForChr[words[0]].append(int(words[2]))
            depth = int(words[2])
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

    for key in outStats.keys():
        if len(coords[key]) > 0:
            outStats[key]['region_start'] = min(coords[key])
            outStats[key]['region_end'] = max(coords[key])
            outStats[key]['region_length'] = outStats[key]['region_end']-outStats[key]['region_start']+1
        if len(depthsForChr[key]) > 0:
            outStats[key]['average_depth'] = round(mean(depthsForChr[key]),2)
            outStats[key]['median_depth'] = median(depthsForChr[key])

    
    with open(out_file_path,"w") as out:
        out.write("CHR"+"\t"+"region_start"+"\t"+"region_end"+"\t"+"region_length"+"\t"+"n_bp_0X"+"\t"+"n_bp_1X"+"\t"+"n_bp_2X"+"\t"+"n_bp_3X"+"\t"+"n_bp_4X"+"\t"
                    +"n_bp_5X"+"\t"+"n_bp_10X"+"\t"+"n_bp_20X"+"\t"+"n_bp_30X"+"\t"+"n_bp_40X"+"\t"+"n_bp_50X"+"\t"+"n_bp_60X"+"\t"+"n_bp_70X"+"\t"+"n_bp_80X"+"\t"
                    +"n_bp_90X"+"\t"+"n_bp_100X"+"\t"+"average_depth"+"\t"+"median_depth"+"\n")

        for key in outStats.keys():
            if (re.match("chr[1-9][0-2]?$",key)):
                out.write(key+"\t"+str(outStats[key]["region_start"])+"\t"+str(outStats[key]["region_end"])+"\t"+str(outStats[key]["region_length"])+"\t"+str(outStats[key]["n_bp_0X"])+"\t"
                    +str(outStats[key]["n_bp_1X"])+"\t"+str(outStats[key]["n_bp_2X"])+"\t"+str(outStats[key]["n_bp_3X"])+"\t"+str(outStats[key]["n_bp_4X"])+"\t"
                    +str(outStats[key]["n_bp_5X"])+"\t"+str(outStats[key]["n_bp_10X"])+"\t"+str(outStats[key]["n_bp_20X"])+"\t"+str(outStats[key]["n_bp_30X"])+"\t"+str(outStats[key]["n_bp_40X"])+"\t"
                    +str(outStats[key]["n_bp_50X"])+"\t"+str(outStats[key]["n_bp_60X"])+"\t"+str(outStats[key]["n_bp_70X"])+"\t"+str(outStats[key]["n_bp_80X"])+"\t"
                    +str(outStats[key]["n_bp_90X"])+"\t"+str(outStats[key]["n_bp_100X"])+"\t"+str(outStats[key]["average_depth"])+"\t"+str(outStats[key]["median_depth"])+"\n")




def main(argv):
    bam_file_path = ''
    bed_file_path = ''
    out_file_path = ''
    all_autosomal = False
    try:
        opts, args = getopt.getopt(argv,"ai:b:o:")
    except getopt.GetoptError:
        sys.exit(2)
    if not opts:
        sys.exit(2)
    for opt, arg in opts:
        if (opt == "-i"):
            bam_file_path = arg
        elif (opt == "-b"):
            bed_file_path = arg
        elif (opt == "-o"):
            out_file_path = arg
        elif (opt == "-a"):
            all_autosomal = True
    if (all_autosomal):
        generateAllDepth(bam_file_path,out_file_path)
    else:
        generateDepth(bam_file_path,bed_file_path,out_file_path)


if __name__ == '__main__':
    main(sys.argv[1:])