# -*- coding: utf-8 -*-
"""
Finds locations of CTCF binding peaks within <threshold> bp of a loop anchor site of interest and returns
a table with the loop name, which side of the loop the anchor of interest is on, and a string with the 
locations of competing CTCFs

Usage: python3 find_competing_ctcfs.py <loopname_file> <ctcf_peak_file> <loop_location_file> <threshold> <output filename>

@author: Kara Quaid
"""

import sys

#this function takes a single loop name, the dictionary of loop locations, the 
#list of ctcf chromosomes, list of the end positions of ctcfs, and the threshold distance,
#and returns a string with the locations of competing ctcf peaks in the format:
#Lx1,Lx2,... where x1 is the distance (bp) from the loop anchor of the first competing
#ctcf, x2 is the distance for the second competing ctcf, etc. 
def searchleft(loopname,loopdict,cchrs,cends,threshold):
    #This function finds ctcfs that are within the threshold distance from the 
    #left end of the loop anchor
    lch = loopdict[loopname][0]
    lstart = loopdict[loopname][1]
    conflicts = "LNone"
    firstconflict = 1
    min_end = lstart - threshold
    for i in range(len(cends)):
        if cchrs[i] == lch:
            if (cends[i] > min_end) & (cends[i] < lstart):
                if firstconflict == 1:
                    conflicts = "L%s"%(lstart-cends[i])
                    firstconflict = 0
                else:
                    conflicts += ",L%s"%(lstart-cends[i])                
    return conflicts

#this function takes a single loop name, the dictionary of loop locations, the 
#list of ctcf chromosomes, list of the start positions of ctcfs, and the threshold distance,
#and returns a string with the locations of competing ctcf peaks in the format:
#Rx1,Rx2,... where x1 is the distance (bp) from the loop anchor of the first competing
#ctcf, x2 is the distance for the second competing ctcf, etc. 
def searchright(loopname,loopdict,cchrs,cstarts,threshold):
    #This function finds ctcfs that are within the threshold distance from the 
    #right end of the loop anchor
    lch = loopdict[loopname][0]
    lend = loopdict[loopname][2]
    conflicts = "RNone"
    firstconflict = 1
    max_start = lend + threshold
    for i in range(len(cends)):
        if cchrs[i] == lch:
            if (cstarts[i] < max_start) & (cstarts[i] > lend):
                if firstconflict == 1:
                    conflicts = "R%s"%(cstarts[i] - lend)
                    firstconflict = 0
                else:
                    conflicts += ",R%s"%(cstarts[i] - lend)               
    return conflicts

##############Main Script##############

loopfile = sys.argv[1]
ctcffile = sys.argv[2]
loop_loc_file = sys.argv[3]
threshold = int(sys.argv[4])
output = sys.argv[5]
#loopfile = "looplist.txt"
#ctcffile = "CTCF_peaks.bed"
#loop_loc_file = "orthologyAnnotated_REderived_human_loop_anchorCTCF.txt"
#threshold = 50000

looplist = []#list of loopnames for loops of interest
with open(loopfile) as f:
    for line in f:
        sline = line.rstrip().split()
        looplist.append(sline[0])

loopdict = {}#dictionary with loopnames as keys and values [chr,start,end]
with open(loop_loc_file) as f:
    for line in f:
        sline = line.rstrip().split()
        loopname = sline[3]
        ch = sline[0]
        start = int(sline[1])
        end = int(sline[2])
        side = sline[5]
        if loopname in looplist:
            newname = loopname + "_" + side
            loopdict[newname] = [ch,start,end,loopname,side]

                

cchrs = []#list of chromosomes for ctcf peaks
cstarts = [] #start positions for ctcf peaks
cends = []#end positions for ctcf peaks
with open(ctcffile) as f:
    firstline = 1
    for line in f:
        if firstline == 1:
            firstline = 0
        else:
            sline = line.rstrip().split()
            cchrs.append(sline[0])
            cstarts.append(int(sline[1]))
            cends.append(int(sline[2]))


conflict_dict = {}
#make a dictionary with loopnames as keys, and strings with the positions of competing
#ctcfs as values
for key in loopdict:
    l_conflicts = searchleft(key,loopdict,cchrs,cends,threshold)
    r_conflicts = searchright(key,loopdict,cchrs,cstarts,threshold)
    all_conflicts = l_conflicts +","+ r_conflicts
    conflict_dict[key] = [all_conflicts,loopdict[key][3],loopdict[key][4]]
    
#create a tab-delimited file with columns: loopname, re-derived_side, competing ctcf string
with open(output,"w") as nf:
    nf.write("loop_name\tre_side\tconflicting_ctcf\n")
    for key in conflict_dict:
        nf.write("%s\t%s\t%s\n"%(conflict_dict[key][1],conflict_dict[key][2],conflict_dict[key][0]))

