from Bio import SeqIO
from Bio.SeqUtils import GC
import numpy as np
import glob

glob_AS = glob.glob('./ncbi_antismash/nf_output/*/')
genome_list,taxa_list,fragBGC,compBGC,totalBGC,seen = [],[],[],[],[],[]

for AS_path in glob_AS:
    glob_BGCs = glob.glob('%s/*region*.gbk'%AS_path)
    if glob_BGCs:
        frag_count = 0
        comp_count = 0
        for BGC in glob_BGCs:
            input_handle = open(BGC,"r")
            for record in SeqIO.parse(input_handle,"genbank"):
                for feat in record.features:
                    if feat.type == 'cand_cluster':
                        if BGC not in seen:
                            seen.append(BGC)
                            if str(feat.qualifiers['contig_edge']) == True:
                                frag_count += 1
                            else:
                                comp_count += 1
        fragBGC.append(frag_count)
        compBGC.append(comp_count)
        totalBGC.append(frag_count+comp_count)
    else:
        fragBGC.append(0)
        compBGC.append(0)
        totalBGC.append(0)
        
print(len(fragBGC),len(compBGC),len(totalBGC))