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


def read_bed(peak_path):
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
            peak[3] = str( peak[3] )
            peak[4] = str( peak[4] )
            peak[5] = str( peak[5] )
            peak = peak[0:6]
            peaks.append(peak)

    return peaks

def check_tasks(task_list):

    possible_tasks = ['DATA', 'BUILD', 'PREDICT']
    task_list.sort()
    
    if ( len(task_list) == 3 ):
        print("can't use three tasks at a time")
        exit(2)

    if ( len(task_list) == 2 ) and \
       ( task_list != ["BUILD", "DATA"] ):
        print("Combination of task not allowed")
        exit(2)

    if( task_list[0] not in possible_tasks ):
        print(task_list[0], "not a TASK(s) in [DATA, BUILD or PREDICT]")
        exit(2)



def check_arguments(arguments):

    if "DATA" in arguments.TASK:

        try:
            int(arguments.RANGE)
        except:
            print("--range should be an integer")
            exit(2)
        try:
            int(arguments.TAG_SUMMIT_RANGE)
        except:
            print("--tag_summit_range should be an integer")
            exit(2)

        if not os.path.isfile(arguments.BAM_CONTROL):
            print("--bam_control file doesn't exist")
            exit(2)

        if not os.path.isfile(arguments.BAM_PULLDOWN):
            print("--bam_pulldown file doesn't exist")
            exit(2)
            
        if not os.path.isfile(arguments.PEAKS):
            print("--regions file doesn't exist")
            exit(2)

        if not os.path.isfile(arguments.PFM):
            print("--pfm file doesn't exist")
            exit(2)

        if not os.path.isfile(arguments.MET):
            print("--wgbs_met_data file doesn't exist")
            exit(2)

        if not os.path.isfile(arguments.UNMET):
            print("--wgbs_unmet_data file doesn't exist")
            exit(2)

        if not os.path.isfile(arguments.DNA_ACC):
            print("--dna_acc_map file doesn't exist")
            exit(2)

        if not os.path.isfile(arguments.GENOME):
            print("--genome_fa file doesn't exist")
            exit(2)
            
        if not os.path.isfile(arguments.CHR_SIZES):
            print("--chr_sizes file doesn't exist")
            exit(2)

        if not os.path.isfile(arguments.REPEATS):
            print("--mask_regions file doesn't exist")
            exit(2)

        if not os.path.isdir(arguments.IN_PATH):
            print("--data_dir directory doesn't exist")
            exit(2)

        if (arguments.OUT_PATH == "default_none"):
            print("--output_dir you should specify an output directory")
            exit(2)
            
            
            
    if "BUILD" in arguments.TASK:

        try:
            int(arguments.FLANKING)
        except:
            print("--flanking should be an integer")
            exit(2)

        if not os.path.isdir(arguments.IN_PATH):
            print("--data_dir directory doesn't exist")
            exit(2)

        if (arguments.OUT_PATH == "default_none"):
            print("--output_dir you should specify an output directory")
            exit(2)


    if "PREDICT" in arguments.TASK:

        if not os.path.isfile(arguments.COEFFS):
            print("--model_coeficients file doesn't exist")
            exit(2)

        if arguments.MOTIF_IN_REGIONS not in ["SEARCH", "START"]:
            print("--motif_in_regions should be either SEARCH or START")
            exit(2)

        if not os.path.isfile(arguments.PEAKS):
            print("--regions file doesn't exist")
            exit(2)

        if not os.path.isfile(arguments.PFM):
            print("--pfm file doesn't exist")
            exit(2)

        if not os.path.isfile(arguments.MET):
            print("--wgbs_met_data file doesn't exist")
            exit(2)

        if not os.path.isfile(arguments.UNMET):
            print("--wgbs_unmet_data file doesn't exist")
            exit(2)

        if not os.path.isfile(arguments.DNA_ACC):
            print("--dna_acc_map file doesn't exist")
            exit(2)

        if not os.path.isfile(arguments.GENOME):
            print("--genome_fa file doesn't exist")
            exit(2)
            
        if not os.path.isfile(arguments.CHR_SIZES):
            print("--chr_sizes file doesn't exist")
            exit(2)

        if not os.path.isdir(arguments.IN_PATH):
            print("--data_dir directory doesn't exist")
            exit(2)
            
        try:
            int(arguments.FLANKING)
        except:
            print("--flanking should be an integer")
            exit(2)

        if (arguments.OUT_PATH == "default_none"):
            print("--output_dir you should specify an output directory")
            exit(2)
            
            
