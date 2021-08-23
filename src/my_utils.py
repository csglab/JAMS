import os
import subprocess

chromosomes = ["chr1", "chr2", "chr3", "chr4", "chr5", \
               "chr6", "chr7", "chr8", "chr9", "chr10", \
               "chr11", "chr12", "chr13", "chr14", "chr15", \
               "chr16", "chr17", "chr18", "chr19", "chr20", \
               "chr21", "chr22", "chrX", "chrY", "chrM"]

def run_cmd(cmd, dry=False):

    if not dry:
        print(cmd)
        # subprocess.call(cmd, shell=True)
        subprocess.run(cmd, shell=True) #, check=True)
        
    else:
        print("## Dry:")
        print(cmd)

def check_path(path):
    # Is path global or relative?
    if path[0] == ".":
        return "%s/%s" % (os.getcwd(), path)
    else:
        return path
    
def remove_tmp_files(out_path):

    rm_cmd = "find %s -type f -name '*.tmp' -delete" % (out_path)
    run_cmd(rm_cmd)
    
def check_arguments(args):
    print("Check formats TODO")
    ## Check if files exist, args are strings, formats are right,
    ##  move this to function script
    return True
    
def read_line(line):
    new_line = line.replace("\n", "").split("\t")
    return new_line


def read_chr_size(chr_path):
    
    with open(chr_path, "r") as sizes:
        chr_sizes = {}
        for chrom in sizes:
            chrom = read_line(chrom)
            chr_sizes[chrom[0]] = int(chrom[1])

    return chr_sizes


def read_peaks(peak_path):
    peaks = []
    
    with open(peak_path, "r") as peak_line: 
        for peak in peak_line:

            if peak[0] == "#": # Skip comments
                continue

            peak = read_line(peak)

            if len(peak) == 1 or peak[0] == "chr":
                continue
            
            peak[0] = str( peak[0] )

            peak[1] = int( peak[1] )
            peak[2] = int( peak[2] )
            peak[3] = int( peak[3] )
            peak[4] = int( peak[4] )
            peak[5] = float( peak[5] )
            peak[6] = float( peak[6] )
            peak[7] = float( peak[7] )
            peak[8] = float( peak[8] )
            peak = peak[0:9]
            # For MACS2 files, 
            # convert abs_summit to summit.
            if peak[4] >= peak[1]:
                peak[4] = peak[4] - peak[1]                

            peaks.append(peak)

    return peaks



def read_encode(encode_path):
    peaks = [ ]
    with open(encode_path, "r") as peak_line: 

        for peak in peak_line:

            if peak[0] == "#":
                continue            

            peak = read_line(peak)
            peaks.append(peak)
            
    return peaks

def out_of_chrom(peak, chr_sizes):

    chr_sizes = read_chr_size(chr_sizes)

    try:
        if int(peak[1]) <= 0 or \
           int(peak[2]) > chr_sizes[peak[0]] or \
           peak[0] not in chromosomes:

            # print("out of chromosome")
            return True
    
        else:
            return False

    except KeyError as e:
        return False

