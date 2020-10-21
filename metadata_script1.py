from Bio import SeqIO
from Bio.SeqUtils import GC
import numpy as np
import glob
import pandas as pd

glob_AS = glob.glob('./ncbi_antismash/nf_output/*/')
genome_list,taxa_list,fragBGC,compBGC,totalBGC,seen = [],[],[],[],[],[]

for AS_path in glob_AS:
    genome_file = (AS_path.split('/')[3]).split('_')[0]
    genome_list.append(genome_file.split('.')[0])
    genome_path = "./ncbi_genomes/%s"%genome_file
    with open(genome_path) as f:
        first_line = f.readline()
        if 'TPA_asm:' in str(first_line):
            taxa = first_line.split(' ')[2]
        if 'uncultured' in str(first_line):
            taxa = first_line.split(' ')[2]
        if 'TPA_asm:' not in str(first_line) and 'uncultured' not in str(first_line):
            taxa = first_line.split(' ')[1]
        taxa_list.append(taxa)
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
                            if str(feat.qualifiers['contig_edge'][0]) == 'True':
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
    

def get_draft_counts(query_genome):
    input_handle = open(query_genome,"r")
    node_count = 0
    gc_count = []
    for record in SeqIO.parse(input_handle,"fasta"):
        node_count += 1
        gc_count.append(GC(record.seq))
    gc_average = np.average(gc_count)
    return node_count,gc_average
    input_handle.close()

node_list,gc_list = [],[]
    
for genome in genome_list:
    node_count,gc_average = get_draft_counts('./ncbi_genomes/%s.fasta'%genome)
    node_list.append(node_count)
    gc_list.append(gc_average)

print(len(genome_list),len(taxa_list),len(node_list),len(gc_list),len(fragBGC),len(compBGC),len(totalBGC))

frames = {'GenomeID':genome_list,'Taxa':taxa_list,'Scaffold_count':node_list,
         'GC_content':gc_list,'Fragmented_BGCs':fragBGC,'Complete_BGCs':compBGC,
         'Total_BGCs':totalBGC}

metadata_df = pd.DataFrame(data=frames)

metadata_df.to_csv('cyanobiome_metadata_df-TFL200507.tsv',sep='\t')