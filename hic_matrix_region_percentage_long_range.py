# -*- coding: utf-8 -*-
"""
Created on Thu Jan 21 13:24:42 2021

@author: 15017
"""
#input dumped matrix file and two regions of interest that the reads will fall between
import sys
filename=sys.argv[1];
region1_start= int(sys.argv[2]);
region1_end=int(sys.argv[3]);
region2_start=int(sys.argv[4]);
region2_end=int(sys.argv[5]);

#enriched_start=157100000
#starting at affected region
enriched_start=158400000
#enriched_end=159100000
#cut off last 100kb
enriched_end=159000000
inside_region=0
total=0
with open(filename) as f:
    for line in f:
        sline = line.strip().split()
        start=int(sline[0]);
        end=int(sline[1]);
        count=float(sline[2]);
        if (end-start)> 30000:
            if (start>=enriched_start) & (end <= enriched_end):
                total=total+count;
                if (start > region1_start) & (start < region1_end):
                    if (end > region2_start) & (end < region2_end):
                        inside_region=inside_region+count
                elif (start > region2_start) & (start < region2_end):
                    if(end > region1_start) & (end < region1_end):
                        inside_region=inside_region+count
                
                
print(total)
print(inside_region)
print(inside_region/total)   
