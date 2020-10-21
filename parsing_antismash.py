import pandas as pd
import os
import subprocess
import re
from Bio import SeqIO

def split_string_at(s, c, n):
    words = s.split(c)
    return c.join(words[:n]), c.join(words[n:])

def download_mibig(bigscape_inputs):
    print("Downloading and parsing MIBiG database of previously characterized BGCs")
    tar_file = os.path.join(bigscape_inputs,"mibig_gbk_2.0.tar.gz")
    tar_cmd = "wget https://dl.secondarymetabolites.org/mibig/mibig_gbk_2.0.tar.gz -O %s"%(tar_file)
    subprocess.call(tar_cmd,shell=True)
    subprocess.call("tar xvzf %s -C %s"%(tar_file,bigscape_inputs),shell=True)
    subprocess.call("rm %s"%tar_file,shell=True)
    
def parse_gbk_list(bigscape_inputs,folder_list):
    gbk_list,new_name_list = [],[]
    for antismash_folder in folder_list:
        print("Parsing BGCs from folder %s"%antismash_folder)
        for root, dirs, files in os.walk(antismash_folder):
            count = 1
            for file in sorted(files):
                if file.endswith(".gbk"):
                    if 'region' in file:
                        strain_name = os.path.basename(os.path.normpath(root)).split('.')[0]
                        if len(strain_name) == 6:
                            gbk_name = file.split('.')[0][:6]
                        else:
                            gbk_name = file.split('.')[0][:8]
                        if gbk_name == strain_name:
                            gbk_list.append(os.path.join(root, file))
                            new_name = file.rstrip('.gbk')
                            new_name_list.append(new_name)
                        else:
                            new_name = strain_name + '.region' + "{0:0=3d}".format(count)
                            count += 1
                            new_name_list.append(new_name)
                            gbk_list.append(os.path.join(root,file))
    return gbk_list,new_name_list

def rename_antismash_gbk(new_name,filepath,outpath):
    new_file = []
    input_handle = open('%s'%filepath, "r")
    for record in SeqIO.parse(input_handle, "genbank"):
        record.id = new_name
        record.description = new_name
        new_file.append(record)
    output_handle = open("%s/%s.gbk"%(outpath,new_name), "w")
    SeqIO.write(new_file, output_handle, "genbank")
    output_handle.close()
    input_handle.close()

def parsing_antismash(bigscape_inputs,folder_list):
    subprocess.call("mkdir -p %s"%(bigscape_inputs),shell=True)
    download_mibig(bigscape_inputs)
    gbk_list,new_name_list = parse_gbk_list(bigscape_inputs,folder_list)
    for i,item in enumerate(gbk_list):
        new_name = new_name_list[i]
        rename_antismash_gbk(new_name,item,bigscape_inputs)
        
bigscape_inputs = "/home/jovyan/jupyterdata/tiago/cyanobiome/bigscape-inputs/"
folder_list = ["/home/jovyan/jupyterdata/tiago/cyanobiome/ncbi_antismash/nf_output/"]

parsing_antismash(bigscape_inputs,folder_list)