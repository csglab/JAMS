

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

            if peak_format == "MACS2":
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
