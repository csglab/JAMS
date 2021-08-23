from . import my_utils
import os

def run_rcade2(tools, args):

    rcade_out = args.OUT_PATH + "/rcade2"
    rcade_prefix = rcade_out + "/" + args.EXP_ID
    my_utils.run_cmd("rm -R %s " % (rcade_out))
    my_utils.run_cmd("mkdir -p %s" % (rcade_out))
    top_peaks = "%s.top_%s_peaks.bed" % (rcade_prefix, args.RC_TOP_PEAKS)

    get_top_peaks(args.PEAKS, args.RC_TOP_PEAKS, \
                  rcade_prefix, args.CHR_SIZES, args.REPEATS, \
                  top_peaks, args.RC_SUMMIT)

    peaks_seq = "%s.top_%s_summits.%sbp.fasta" % \
                (rcade_prefix, str(args.RC_TOP_PEAKS), str(args.RC_SUMMIT))
    
    ## Get peak protein sequences
    peak_seq_cmd = "bedtools getfasta -fi %s -bed %s > %s" % \
                   (args.GENOME, top_peaks, peaks_seq)
    my_utils.run_cmd(peak_seq_cmd)

    ## check if path is global
    c2h2_seq = my_utils.check_path(args.C2H2_SEQ)
    peaks_seq = my_utils.check_path(peaks_seq)
    
    # Run RCADE2
    rcade_cmd = "cd %s; bash RCOpt.sh %s %s %s" % (tools["RCADE2"], \
                                                   args.EXP_ID, \
                                                   c2h2_seq, peaks_seq)
    my_utils.run_cmd(rcade_cmd)
    # Here is where rcade dumps the results.
    rcade_path = "%s/out/%s" % (tools["RCADE2"], args.EXP_ID)
    mv_cmd = "mv %s %s" %(rcade_path, rcade_out)
    my_utils.run_cmd(mv_cmd)


    PFM = "%s/%s/results.opt.PFM.txt" %(rcade_out, args.EXP_ID)
    
    # Check if exists, 
    if not os.path.isfile(PFM):
        print("Could not create the PFM for C2H2 protein. Check RCADE log")
    
    return PFM


def get_top_peaks(peaks, top, rcade_prefix, chr_sizes, repeats, top_out, \
                  rc_summit_region):

    ## Format peak bed, eliminate the out of chromosomes entries.
    #   get peak at +- 250 region from the start of the summit.
    peaks_tmp = rcade_prefix + ".tmp"
    with open(peaks_tmp, "w") as peaks_tmp_f:

        for peak in peaks:

            summit = peak[4]
                
            new_start = str(int(peak[1]) + int(summit) - \
                            int(float(rc_summit_region)/2))

            new_stop =  str(int(peak[1]) + int(summit) + \
                            int(float(rc_summit_region)/2))
            
            new_line = [peak[0], new_start, new_stop, ".", 
                        str( peak[6] ) ]

            if my_utils.out_of_chrom(new_line, chr_sizes):
                continue

            # print(new_line)
            new_line = "\t".join(new_line)
            
            print(new_line, sep="\t", end="\n", file=peaks_tmp_f,
                  flush=False)

    ## Sort and then mask repetitive sequences.
    sort_mask_cmd ="sort -k1,1 -k2,2n %s | bedtools intersect -sorted -v -a - -b %s | sort -k 5,5nr  - | head -n %s > %s" \
      % (peaks_tmp, repeats, top, top_out)

    my_utils.run_cmd(sort_mask_cmd)
    
