################################################
##############   FUNCTION    ###################
################################################
from . import my_utils

def write_summit_and_score_beds(summit_path, score_path,  \
                                peaks, out_path, \
                                chr_sizes, repeats, \
                                this_range, peak_format):

    tmp_path = out_path + "/bed_without_repeats.bed.tmp"
    with open(tmp_path, "w") as tmp_out:
        
        for peak in peaks:
            
            chrom = peak[0]
            old_start_peak = peak[1]
            old_stop_peak = peak[2]

            name_old_peaks = "%s:%s-%s" % (chrom, old_start_peak, \
                                           old_stop_peak)

            if peak_format == "MACS":
                new_start = int(old_start_peak) + int(peak[4]) - \
                            int(this_range)
                new_stop = int(old_start_peak) + int(peak[4]) + \
                           int(this_range)

                tag = peak[5]
                length = peak[3]
                score = peak[6]
                enrichment = peak[7]
                
            if peak_format == "ENCODE":
                new_start = int(old_start_peak) + int(peak[9]) - \
                            int(this_range)
                new_stop = int(old_start_peak) + int(peak[9]) + \
                           int(this_range)
                tag = peak[10] # peak[-1]
                length = int(old_stop_peak) - int(old_start_peak) + 1
                score = peak[8] # peak[8], using pval7, is qval +?
                enrichment = peak[6]

            name_new_peaks = "%s:%s-%s" % (chrom, new_start, \
                                           new_stop)
            
            name = "%s::%s" % (name_old_peaks, name_new_peaks)

            new_line = [chrom, str(new_start), str(new_stop), \
                        name, str(tag), str(length), str(score), \
                        str(enrichment), name_old_peaks]
            
            if my_utils.out_of_chrom(new_line, chr_sizes):
                continue
            new_line = "\t".join(new_line)
            print(new_line, end="\n", file=tmp_out)

    tmp_path2 = tmp_path + ".tmp"
    bed_cmd = "bedtools intersect -v -a %s -b %s > %s" \
              % (tmp_path, repeats, tmp_path2)

    my_utils.run_cmd(bed_cmd)
    
    with open(tmp_path2, "r") as tmp_file, \
         open(summit_path, "w") as summit, \
         open(score_path, "w") as score:

        # Write header to summit file
        header = "Name\tfold_enrichment\ttags\tlength\tMACS_score"
        print(header, end="\n", file=score)
        
        for peak in tmp_file:
            
            peak = my_utils.read_line(peak)
            # chr1    15521   15721   chr1:15120-15904        5.51
            summit_line = "%s\t%s\t%s\t%s\t%s" % (peak[0],
                                                  peak[1],
                                                  peak[2],
                                                  peak[-1],\
                                                  peak[6])
            print(summit_line, end="\n", file=summit)
            
            # Name    fold_enrichment tags    length  MACS_score
            # chr1:15120-15904::chr1:15521-15721  5.51  39  785  63.7
            score_line = "%s\t%s\t%s\t%s\t%s" % \
                         (peak[3], \
                          peak[7], \
                          peak[4], \
                          peak[5], \
                          peak[6])
            print(score_line, end="\n", file=score)

def get_fasta_file(genome, bed_summits, fasta_file):

    tmp_bed_file = bed_summits + ".tmp"
    ## Create tmp bed file, with name: third column::chr:start-stop
    with open(bed_summits, "r") as original_bed, \
         open(tmp_bed_file, "w") as tmp_bed:

        for peak in original_bed:
            peak = my_utils.read_line(peak)

            new_name = "%s::%s:%s-%s" \
                       % (peak[3], peak[0], peak[1], peak[2])
            new_line = [peak[0], peak[1], peak[2], new_name, peak[4]]
            new_line = "\t".join(new_line)
            print(new_line, end="\n", file=tmp_bed)
    
    get_fasta_cmd = "bedtools getfasta -name -fi %s -bed %s -fo %s" % \
                    (genome, tmp_bed_file, fasta_file)

    my_utils.run_cmd(get_fasta_cmd)

def write_aligned_positions_bed(affimx_position_path, \
                                aligned_positions_bed):
    
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

    # bed_cmd = "bedtools getfasta -name -s -tab -fi %s -bed %s -fo %s" % \
    #           (genome, tmp_out, tmp_out2)
    bed_cmd = '''bedtools getfasta -name -s -tab -fi %s -bed %s | awk 'BEGIN {{OFS = "\\t"}} $2 = toupper($2)' - > %s''' % (genome, tmp_out, tmp_out2)
    
    my_utils.run_cmd(bed_cmd)

    with open(tmp_out, "r") as tmp_f, \
         open(tmp_out2, "r") as tmp2_f, \
         open(aligned_fasta, "w") as fasta:

        for line1, line2 in zip(tmp_f, tmp2_f):

            line1, line2 = my_utils.read_line(line1), my_utils.read_line(line2)

            if not line1[3] == line2[0][:-3]:
                print("Error")
                exit

            new_line = [line1[3], line1[0], line1[1],\
                        line1[2],  line1[5], line2[1]]

            new_line = "\t".join(new_line)
            print(new_line, end="\n", file=fasta)

def write_aligned_numeric_tab(aligned_fasta, aligned_numeric, aligned_tab):

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


def get_methylreads(aligned_positions_bed, \
                    methyl_map, \
                    seq_methyl):

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
