from . import my_utils
import os

def extract_data(args):
    """
    """
    MC_INPUT_DIR = args.IN_PATH
    # define prefix for files
    out_prefix = "%s/%s" % \
        (MC_INPUT_DIR, args.EXP_ID)
    
    # Remove dir from previous run.
    # my_utils.run_cmd("rm -R %s" % (MC_INPUT_DIR))
    my_utils.run_cmd("mkdir -p %s" % (MC_INPUT_DIR))

    # Copy PFM to input directory.
    if args.PFM != "None":
        my_utils.run_cmd( "cp %s %s" % \
                          (args.PFM, MC_INPUT_DIR))

    ## 01. Create summit files
    if "BUILD" in args.TASK:
    
        write_summit_and_score_beds(out_prefix, args.PEAKS, \
                                    MC_INPUT_DIR, args.CHR_SIZES, \
                                    args.REPEATS, args.RANGE, args)
    elif "PREDICT" in args.TASK:
        
        ## Get motif length
        pfm_length = int(get_PFM_length_from_model(args.COEFFS))

        if "SEARCH" == args.MOTIF_IN_REGIONS:
            ## write bed and summit vRepeat files with
            ##  the same coordinates
            write_summit_and_score_beds_predt(out_prefix, \
                                              args.PEAKS, \
                                              args.CHR_SIZES)
    
    ##### 03.Extract sequence for summits
    if "START" != args.MOTIF_IN_REGIONS:
        print("Retrieving summit sequences ... ")
        get_fasta_file(args.GENOME, out_prefix)

    aligned_positions_bed = "%s_aligned_positions.bed" % (out_prefix)
    
    if ( "BUILD" in args.TASK ) or \
       ( "PREDICT" in args.TASK and \
         "SEARCH" == args.MOTIF_IN_REGIONS ):
        
        ##### 04. Run Affimx
        print("Running AffiMx ... ")
        # From task 3
        fasta_path= "%s_summits_vRepeats.fa" % (out_prefix)
        affimx_prefix = out_prefix + "_affimx"
        affimx_cmd = "%s/AffiMx -pwm %s -fasta %s -out %s > %s.log" % \
            (os.path.dirname(os.path.realpath(__file__)), \
             args.PFM, fasta_path, affimx_prefix, affimx_prefix)

        my_utils.run_cmd(affimx_cmd)

        ##### 05. write_aligned_positions_bed
        # Output from Affimx
        affimx_position_path = "%s.position.txt" % (affimx_prefix)
                
        #  Create the BED file of aligned genomic coordinates
        write_aligned_positions_bed(affimx_position_path, \
                                    aligned_positions_bed)

    elif( "PREDICT" in args.TASK and \
          "START" == args.MOTIF_IN_REGIONS ):
        write_aligned_pos_bed_n_txt(out_prefix, args.PEAKS, args.CHR_SIZES)
        ## write aligned file
        pass
    
    
    ##### 06.write aligned fasta
    aligned_fasta = out_prefix + "_aligned_sequences.txt"
    write_aligned_fasta(aligned_positions_bed, \
                        aligned_fasta, args.GENOME, \
                        args.CHR_SIZES, args.RANGE)
    
    ##### 07.write aligned numeric tab
    aligned_numeric = out_prefix + "_aligned_sequences_numeric_mx.txt"
    aligned_tab = out_prefix + "_aligned_sequences_tabulated_mx.txt"
    write_aligned_numeric_tab(aligned_fasta, \
                              aligned_numeric, \
                              aligned_tab)
    ##### 08.
    ### Nonmethylated 
    nonmet_seq = out_prefix + "_aligned_sequences.fasta.nonmethylreads"
    get_methylreads(aligned_positions_bed, args.UNMET, nonmet_seq)
    ##### 09.
    ### Methylated
    met_seq = out_prefix + "_aligned_sequences.fasta.methylreads"
    get_methylreads(aligned_positions_bed, args.MET, met_seq)
    ##### 10.
    ## DNase-seq
    dna_acc = out_prefix + "_aligned_sequences.fasta.accessibility"
    get_dna_acc(aligned_positions_bed, args.DNA_ACC, dna_acc)

    
    if "BUILD" in args.TASK:
        ##### 10. Extract tags from control and pulldown bams
        extract_tags_ctrl_pd(out_prefix,
                             args.BAM_CONTROL,
                             args.BAM_PULLDOWN,
                             args.TAG_SUMMIT_RANGE)

    my_utils.remove_tmp_files(args.IN_PATH)

    
def run_methyl_chip(args):

    MC_INPUT_DIR = args.IN_PATH
    methyl_dir = my_utils.check_path(MC_INPUT_DIR)
    my_utils.remove_tmp_files(MC_INPUT_DIR)

    methyl_cmd = "Rscript %s/JAMS_GLM.R --experiment %s --flanking %s --input_dir %s --output_dir %s --script_path %s" % \
        ( os.path.dirname(os.path.realpath(__file__)), \
          args.EXP_ID, \
          args.FLANKING, \
          MC_INPUT_DIR, \
          args.OUT_PATH, \
          os.path.dirname(os.path.realpath(__file__)) )
        
    my_utils.run_cmd(methyl_cmd, dry=False)


def write_summit_and_score_beds(out_prefix, peaks, out_path,
                                chr_sizes, repeats, this_range, args):
    """
    Task 1
    """
    # Input: peaks
    # Temporary files.
    formated_peaks_tmp = "%s_T1_F1_formated_summit_peaks.tmp" % (out_prefix)
    no_repeats_bed_tmp = "%s_T1_F2_formated_summit_peaks_no_repeats.tmp" % \
        (out_prefix)
    
    # Output files.
    summit_path = "%s_summits_vRepeats.bed" % (out_prefix)
    score_path = "%s_summits_vRepeats_scores.txt" % (out_prefix)
    
    with open(formated_peaks_tmp, "w") as formated_peaks_f:
        
        for peak in peaks:
            
            chrom = peak[0]
            old_start = peak[1]
            old_stop = peak[2]
            
            ## New start and stop of peak +/- range from summit
            new_start = int(old_start) + int(peak[4]) - int(this_range)
            new_stop = int(old_start) + int(peak[4]) + int(this_range)

            tag = peak[5]
            length = peak[3]
            score = peak[6]
            enrichment = peak[7]
            
            name_old_peaks = "%s:%s-%s" % (chrom, old_start, old_stop)
            name_new_peaks = "%s:%s-%s" % (chrom, new_start, new_stop)
            name = "%s::%s" % (name_old_peaks, name_new_peaks)

            
            new_line = [str(chrom), str(new_start), str(new_stop),
                        name, str(tag), str(length), str(score),
                        str(enrichment), name_old_peaks]

            # Check if peak stop is outside chromosome
            if my_utils.out_of_chrom(new_line, chr_sizes):
                continue

            new_line = "\t".join(new_line)
            print(new_line, end="\n", file=formated_peaks_f)

    ## Mask peaks over repeat (or blacklisted) regions
    mask_bed_cmd = "bedtools intersect -v -a %s -b %s > %s" \
        % (formated_peaks_tmp, repeats, no_repeats_bed_tmp)
    my_utils.run_cmd(mask_bed_cmd)

    
    with open(no_repeats_bed_tmp, "r") as no_repeats_bed_f, \
         open(summit_path, "w") as summit_f, \
         open(score_path, "w") as score_f:

        # Write header to summit file
        header = "Name\tfold_enrichment\ttags\tlength\tMACS_score"
        print(header, end="\n", file=score_f)

        for peak in no_repeats_bed_f:

            ## Write summit file "_summits_vRepeats.bed"
            peak = my_utils.read_line(peak)
            # chr1    15521   15721   chr1:15120-15904        5.51
            summit_line = "%s\t%s\t%s\t%s\t%s" % \
              (peak[0], peak[1], peak[2], peak[-1], peak[6])

            print(summit_line, end="\n", file=summit_f)

            ## Write score file "_summits_vRepeats_scores.txt.tmp"
            # Name    fold_enrichment tags    length  MACS_score
            # chr1:15120-15904::chr1:15521-15721  5.51  39  785  63.7
            score_line = "%s\t%s\t%s\t%s\t%s" % \
                         (peak[3], peak[7], peak[4], peak[5], peak[6])
            
            print(score_line, end="\n", file=score_f)


def get_fasta_file(genome, out_prefix):
    """
    Task 3
    """
    ## Input # From Task 1
    summit_path = "%s_summits_vRepeats.bed" % (out_prefix)
    # Temporary file
    tmp_bed = "%s_T3_F1_vRepeats_bed.tmp" % (out_prefix)
    ## Out file
    fasta_path= "%s_summits_vRepeats.fa" % (out_prefix)
    
    ## Create tmp bed file, with name: third column::chr:start-stop
    with open(summit_path, "r") as summit_f, \
         open(tmp_bed, "w") as tmp_bed_f:

        for peak in summit_f:
            peak = my_utils.read_line(peak)

            new_name = "%s::%s:%s-%s" \
                       % (peak[3], peak[0], peak[1], peak[2])
            new_line = [peak[0], peak[1], peak[2], new_name]
            new_line = "\t".join(new_line)
            print(new_line, end="\n", file=tmp_bed_f)
    
    get_fasta_cmd = "bedtools getfasta -name -fi %s -bed %s -fo %s" % \
                    (genome, tmp_bed, fasta_path)

    my_utils.run_cmd(get_fasta_cmd)

    
def write_aligned_positions_bed(affimx_position_path, \
                                aligned_positions_bed):
    """
    Task 05 
    """
    with open(affimx_position_path, "r") as affimx_pos, \
         open(aligned_positions_bed, "w") as aligned_pos:

        # Skip header
        print(affimx_position_path)
        next(affimx_pos)
        
        for peak in affimx_pos:

            peak = my_utils.read_line(peak)
            old_start = int (peak[0].split("::")[1].split(":")[1]\
                             .split("-")[0])
            chrom = peak[0].split("::")[1].split(":")[0]
            am_num = int(peak[1])
            
            if am_num >= 0:
                direction = "+" 
            else:
                direction = "-"

            new_stop = old_start + abs(am_num)
            new_start = new_stop - 1
            
            new_line = "%s\t%s\t%s\t%s\t0\t%s" % (chrom, new_start, \
                                                  new_stop, \
                                                  peak[0], direction)


            print(new_line, end="\n", file=aligned_pos)

            
def write_aligned_fasta(aligned_positions_bed, \
                        aligned_fasta, \
                        genome, chr_sizes, \
                        this_range):
    """
    Task 06 
    """
    tmp_out = aligned_positions_bed + ".tmp"
    tmp_out2 = tmp_out + ".tmp"

    with open(aligned_positions_bed, "r") as positions_bed, \
         open(tmp_out, "w") as tmp_file:

        for peak in positions_bed:

            peak = my_utils.read_line(peak)
            
            if my_utils.out_of_chrom(peak, chr_sizes):
                continue

            peak[1] = str(int(peak[1]) - int(this_range))
            peak[2] = str(int(peak[2]) + int(this_range))

            # Delete entries that start outside the Chromosome
            if int(peak[1]) <= 0:
                continue
            
            peak = "\t".join(peak)
            print(peak, end="\n", file=tmp_file)

    bed_cmd = '''bedtools getfasta -name -s -tab -fi %s -bed %s | awk 'BEGIN {{OFS = "\\t"}} $2 = toupper($2)' - > %s''' % \
      (genome, tmp_out, tmp_out2)
    
    my_utils.run_cmd(bed_cmd)

    with open(tmp_out, "r") as tmp_f, \
         open(tmp_out2, "r") as tmp2_f, \
         open(aligned_fasta, "w") as fasta:

        for line1, line2 in zip(tmp_f, tmp2_f):

            line1, line2 = my_utils.read_line(line1), my_utils.read_line(line2)

            if line1[3] != line2[0][:-3]:
                print("flag")
                print(line1[3])
                print(line2[0])
                print(line2[0][:-3])
                print("Error")
                print(line1)
                print(line2)
                exit(2)

            new_line = [line1[3], line1[0], line1[1],\
                        line1[2],  line1[5], line2[1]]

            new_line = "\t".join(new_line)
            print(new_line, end="\n", file=fasta)

def write_aligned_numeric_tab(aligned_fasta, aligned_numeric, aligned_tab):
    """
    Task 07 write aligned numeric tab
    """
    with open(aligned_fasta, "r") as fasta, \
         open(aligned_numeric, "w") as numeric_f, \
        open(aligned_tab, "w") as tab_f:

        for peak in fasta:

            peak = my_utils.read_line(peak)

            ## Tabular
            seq_tab = "\t".join(list(peak[-1]))
            tab_line = [peak[0], peak[1], peak[2], peak[3], peak[4], \
                        seq_tab]
            tab_line = "\t".join(tab_line)

            print(tab_line, end="\n", file=tab_f)
            
            ## Numeric
            numeric = seq_tab.replace("A", "0").replace("C", "1").\
                      replace("G", "2").replace("T", "3").replace("N", "-1")

            num_line = [peak[0], peak[1], peak[2], peak[3], peak[4], \
                        numeric]
            num_line = "\t".join(num_line)

            print(num_line, end="\n", file=numeric_f)


def get_methylreads(aligned_positions_bed, methyl_map, seq_methyl):
    """
    Task 08 and task 09 get methyl/nonmethyl reads
    """
    seq_methyl_tmp = seq_methyl + ".tmp"
    bw_cmd = "bwtool matrix 100:101 %s %s %s" % \
             (aligned_positions_bed, \
              methyl_map, seq_methyl_tmp)

    my_utils.run_cmd(bw_cmd)

    with open(aligned_positions_bed, "r") as positions, \
         open(seq_methyl_tmp, "r") as reads, \
         open(seq_methyl, "w") as seq_methyl_f:

        for line_pos, line_read in zip(positions, reads):

            line_pos = my_utils.read_line(line_pos)
            line_read = line_read.replace("\n", "")
            
            new_line = [line_pos[3], line_pos[0], line_pos[1], \
                        line_pos[2], line_pos[5], line_read]

            new_line = "\t".join(new_line)
            print(new_line, end="\n", file=seq_methyl_f)

def get_dna_acc(aligned_positions_bed, dna_acc_map, dna_acc_out):
    """
    Task 10 get dna acc data
    """
    print("extract dna accessibility data from ref")
    # bwtool matrix 100:101 ${bed} ${bw} ${out}
    dna_acc_out_tmp = dna_acc_out + ".tmp"

    bw_cmd = "bwtool matrix 1500:1501 %s %s %s" % \
             (aligned_positions_bed, \
              dna_acc_map, dna_acc_out_tmp)

    my_utils.run_cmd(bw_cmd)

    with open(aligned_positions_bed, "r") as positions, \
         open(dna_acc_out_tmp, "r") as reads, \
         open(dna_acc_out, "w") as seq_acc_f:

        for line_pos, line_read in zip(positions, reads):

            line_pos = my_utils.read_line(line_pos)
            line_read = line_read.replace("\n", "")
            
            new_line = [line_pos[3], line_pos[0], line_pos[1], \
                        line_pos[2], line_pos[5], line_read]

            new_line = "\t".join(new_line)
            print(new_line, end="\n", file=seq_acc_f)

    print("Done: extract dna accessibility data from ref")



def extract_tags_ctrl_pd(out_prefix, bam_control,
                         bam_pulldown, tag_summit_range):
    """
    Task 10: get tags.
    """
    ## Input: Should match the one in task 5
    aligned_bed = "%s_aligned_positions.bed" % (out_prefix)

    ## Temporary files
    aligned_bed_mod_tmp = "%s_T10_F1_range_mod_aligned_positions.bed.tmp" % \
      (out_prefix)
    aligned_bed_tag_tmp = "%s_T10_F2_range_mod_wtags_aligned_positions.bed.tmp" % \
      (out_prefix)
    
    ## Output modify
    score_out_tmp = "%s_summits_vRepeats_scores.txt.tmp" % (out_prefix)
    score_out = "%s_summits_vRepeats_scores.txt" % (out_prefix)
    
    # Redefine bed file with start n stop +/- tag_summit_range from summit.
    with open(aligned_bed, "r") as aligned_bed_f, \
        open(aligned_bed_mod_tmp, "w") as aligned_bed_mod_tmp_f:

        # next(aligned_bed_f) # Ignore header

        for aligned_bed_line in aligned_bed_f:

            aligned_bed_line = my_utils.read_line(aligned_bed_line)

            chrom = str(aligned_bed_line[0])
            new_start = int(aligned_bed_line[1]) - int(tag_summit_range)
            new_stop = int(aligned_bed_line[2]) + int(tag_summit_range)

            
            ## Check if new start is negative
            if new_start < 0:
                continue
            
            print(chrom, str(new_start), str(new_stop),
                  file=aligned_bed_mod_tmp_f, sep="\t")
        
    get_tags_cmd = "bedtools multicov -q 30 -bed %s -bams %s %s > %s" % \
      (aligned_bed_mod_tmp, bam_control, bam_pulldown, aligned_bed_tag_tmp)
    
    my_utils.run_cmd(get_tags_cmd)
    cp_cmd = "cp %s %s" % (score_out, score_out_tmp)
    my_utils.run_cmd(cp_cmd)
    
    
    with open(aligned_bed_tag_tmp, "r") as aligned_bed_tag_tmp_f, \
    open(score_out_tmp, "r") as score_out_tmp_f, \
    open(score_out, "w") as score_out_f:
        
        next(score_out_tmp_f)
        score_out_f.write("Name\tfold_enrichment\ttags\tlength\tMACS_score\tctrl.tag.old\tpulldown.tag.old\n")
        
        for peak_tags_line, peak_score_line in \
          zip(aligned_bed_tag_tmp_f, score_out_tmp_f):
            
            peak_score_line = my_utils.read_line(peak_score_line)
            peak_tags_line = my_utils.read_line(peak_tags_line)
            
            print(str(peak_score_line[0]), str(peak_score_line[1]),
                  str(peak_score_line[2]), str(peak_score_line[3]),
                  str(peak_score_line[4]), str(peak_tags_line[-2]),
                  str(peak_tags_line[-1]),
                  file=score_out_f, sep="\t")


def get_PFM_length_from_model(coeff_file):
    """
    """
    with open(coeff_file, 'r') as file: 
        coeff_string = file.read().replace('\n', '')
        occurrences = coeff_string.count("Motif_CpG_state")
    return(occurrences)

def write_summit_and_score_beds_predt(out_prefix, peaks, chr_sizes):
    """
    """
    # Output files.
    summit_path = "%s_summits_vRepeats.bed" % (out_prefix)
    score_path = "%s_summits_vRepeats_scores.txt" % (out_prefix)

    with open(summit_path, "w") as summit_f, \
         open(score_path, "w") as score_f:

        # Write header to summit file
        header = "Name\tscore\tstrand"
        print(header, end="\n", file=score_f)

        
        for peak in peaks:
            # print(peak)
            chrom = peak[0]
            start = peak[1]
            stop = peak[2]
            name = peak[3]
            score = peak[4]
            strand = peak[5]
            
            new_name = "%s:%s-%s" % (chrom, start, stop)
            both_names = "%s::%s" % (new_name, new_name)
            
            ## vRepeats.bed
            summit_line = [str(chrom), str(start), str(stop),
                           str(new_name), str(score), str(strand)]

            score_line = [str(both_names), str(score), str(strand) ]


            # Check if peak stop is outside chromosome
            if my_utils.out_of_chrom(summit_line, chr_sizes):
                continue

            ## Write 
            summit_line = "\t".join(summit_line)
            print(summit_line, end="\n", file=summit_f)

            score_line = "\t".join(score_line)
            print(score_line, end="\n", file=score_f)

            
def write_aligned_pos_bed_n_txt(out_prefix, peaks, chr_sizes):
    """
    """
    # Output files.
    score_path = "%s_summits_vRepeats_scores.txt" % (out_prefix)
    affimx_path = "%s.affimx.affinity.txt" % (out_prefix)
    aligned_path = "%s_aligned_positions.bed" % (out_prefix)

    with open(score_path, "w") as score_f, \
         open(affimx_path, "w") as affimx_f, \
         open(aligned_path, "w") as aligned_f:

        print("Name\tscore\tstrand", end="\n", file=score_f)
        print("Name\tPlaceholder", end="\n", file=affimx_f)

        for peak in peaks:
            # print(peak)
            chrom = peak[0]
            start = peak[1]
            stop = peak[2]
            name = peak[3]
            score = peak[4]
            strand = peak[5]

            old_name = "%s:%s-%s" % (chrom, start, stop)
            new_name = "%s:%s-%s" % (chrom, start, start+1)
            both_names = "%s::%s" % (new_name, old_name)
            
            ## Fake coordinate to check if site is not close to end of chr
            test_line = [str(chrom), str(start), str(start+50)]
            ## Scores
            score_line = [str(both_names), str(score), str(strand) ]
            ## aligned
            aligned_line = [str(chrom), str(start), str(start+1), str(both_names),\
                            "0", str(strand)]
            ## Affimx
            affimx_line = [str(both_names), str("NA") ]

            # Check if peak stop is outside chromosome
            if my_utils.out_of_chrom(test_line, chr_sizes):
                continue

            score_line = "\t".join(score_line)
            print(score_line, end="\n", file=score_f)

            aligned_line = "\t".join(aligned_line)
            print(aligned_line, end="\n", file=aligned_f)

            affimx_line = "\t".join(affimx_line)
            print(affimx_line, end="\n", file=affimx_f)
        
        pass
    
    
