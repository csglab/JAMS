from . import my_utils

def get_top_peaks(peaks, top, rcade_prefix, chr_sizes, repeats, top_out, \
                  rc_summit_region, peak_format):

    ## Format peak bed, eliminate the out of chromosomes entries.
    #   get peak at +- 250 region from the start of the summit.
    peaks_tmp = rcade_prefix + ".tmp"
    with open(peaks, "r") as peaks_f, \
         open(peaks_tmp, "w") as peaks_tmp_f:

        for peak in peaks_f:

            if peak[0] == "#" or peak[0] == "":
                continue
            peak = my_utils.read_line(peak)
            if peak[0] == "chr" or peak == [""]:
                continue
         
            if peak_format == "MACS":
                summit = peak[4]
            if peak_format == "ENCODE":
                summit = peak[9]
                
            new_start = str(int(peak[1]) + int(summit) - \
                            int(float(rc_summit_region)/2))

            new_stop =  str(int(peak[1]) + int(summit) + \
                            int(float(rc_summit_region)/2))
            
            if peak_format == "MACS":
                new_line = [peak[0], new_start, new_stop, ".", \
                            peak[6]]
            
            if peak_format == "ENCODE":
                # qValue (column 9) and then signalValue (column 7)
                new_line = [peak[0], new_start, new_stop, ".", \
                            peak[8], peak[6]]

            if my_utils.out_of_chrom(new_line, chr_sizes):
                continue
            new_line = "\t".join(new_line)
            
            print(new_line, sep="\t", end="\n", file=peaks_tmp_f, flush=False)

    if peak_format == "MACS":
        ## Sort and then mask repetitive sequences.
        sort_mask_cmd ="sort -k1,1 -k2,2n %s | bedtools intersect -sorted -v -a - -b %s | sort -k 5,5nr  - | head -n %s > %s"\
                        % (peaks_tmp, repeats, top, top_out)

    if peak_format == "ENCODE":
        ## Sort and then mask repetitive sequences.
        sort_mask_cmd ="sort -k1,1 -k2,2n %s | bedtools intersect -sorted -v -a - -b %s | sort --reverse --numeric-sort --key=5,5 --key=6,6 - | head -n %s > %s" % (peaks_tmp, repeats, top, top_out)

    my_utils.run_cmd(sort_mask_cmd)
    
