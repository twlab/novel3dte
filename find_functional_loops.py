# -*- coding: utf-8 -*-
"""
Created on Tue Jan  7 10:05:30 2020

@author: Kara
"""

import sys

#reads a bed file of promoters into a dictionary.
def readbed(filename):
    i=0
    dictname = {}
    with open(filename) as bedfile:
        for line in bedfile:
            sline = line.rstrip().split()
            chrom = sline[0]
            start=int(sline[1])
            end=int(sline[2])
            dictname[i]=[chrom,start,end]
            i+=1
    return dictname
        
def find_cre(chromosome,position,cre_dict,threshold,side,bidirectional):#mess around with minval and maxval
    foundstuff={}
    if side=="right":
        maxval=position+threshold
        if bidirectional ==True:
            minval=position - threshold
        else:
            minval=position
    else:
        if bidirectional == True:
            maxval=position + threshold
        else: 
            maxval=position
        minval=position-threshold
    foundone=False
    i=0
    for key in cre_dict:
        cre=cre_dict[key]
        cre_chrom=cre[0]
        if cre_chrom==chromosome:
            cre_start=cre[1]
            if cre_start >= minval:
                cre_end = cre_dict[key][2]
                #print("%s\t%s\t%s"%(cre_chrom,cre_start,cre_end))
                if cre_start<=maxval:
                    foundone=True
                    foundstuff[i]=[cre_chrom,cre_start,cre_end]
                    i+=1
    return foundone,foundstuff
                    
    
#################MAIN SCRIPT###################################################

       
promoter_file="gm_all_promoters.bed"
enhancer_file="gm_polycomb_repressed.bed"
promoters = readbed(promoter_file)
enhancers=readbed(enhancer_file)
loopfile = sys.argv[1]
loop_loc_file = sys.argv[2]
threshold = int(sys.argv[3])
loop_loc_file = "human_loops_re-derived_orthologous.bed"
threshold = 35000
bidirectional = True
#looplist = []#list of loopnames for loops of interest
#with open(loopfile) as f:
#    for line in f:
#        sline = line.rstrip().split()
#        looplist.append(sline[0])
loopdict={}
with open(loop_loc_file) as f:
    for line in f:
        sline = line.rstrip().split()
        loopname = sline[3]
        ch = sline[0]
        start = int(sline[1])
        end = int(sline[2])
        loopdict[loopname] = [ch,start,end,loopname]

candidates_dict={}
cres_dict={}
looplengths=[]
for key in loopdict:
    right_promoters={}
    right_enhancers={}
    left_promoters={}
    left_enhancers={}
    left_anchor=int(loopdict[key][1])
    right_anchor=int(loopdict[key][2])
    chromosome=loopdict[key][0]
    loop_length=right_anchor-left_anchor
    looplengths.append(loop_length)
        #check for promoters within 10kb of right anchor
    right_promoter,right_promoters=find_cre(chromosome,right_anchor,promoters,threshold,"right",bidirectional)
    if right_promoter==True:
        #if found, check for enhancers within 10kb of left anchor
         left_enhancer,left_enhancers=find_cre(chromosome,left_anchor,enhancers,threshold,"left",bidirectional)
         if left_enhancer==True:
            #if found, add loop to dictionary of candidates
            candidates_dict[key]=loopdict[key]
    #check for promoters within 10kb of left anchor
    left_promoter,left_promoters=find_cre(chromosome,left_anchor,promoters,threshold,"left",bidirectional)
    if left_promoter==True:
        #if found, check for enhancers within 10kb of right anchor
        right_enhancer,right_enhancers=find_cre(chromosome,right_anchor,enhancers,threshold,"right",bidirectional)
        if right_enhancer==True:
            #if found, add loop to dictionary of candidates
            if key not in candidates_dict.keys():
                candidates_dict[key]=loopdict[key]
    cres_dict[key]=[right_promoters,left_enhancers,left_promoters,right_enhancers]
                    

with open("human_all_enhancer_promoter_candidates_10kb.txt","w") as nf:
    nf.write("loop_name\tchromosome\tstart\tend\t\n")
    for key in candidates_dict:
        nf.write("%s\t%s\t%s\t%s\t\n"%(key,candidates_dict[key][0],candidates_dict[key][1],candidates_dict[key][2]))

































