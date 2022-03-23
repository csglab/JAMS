#!/bin/bash


Rscript de_novo_motif_discovery_JAMS.R \
	--experiment CTCF_HEK293_de_novo_terminal \
	--flanking 15 \
	--pfm_length 9 \
	--input_dir ./data/CTCF_demo/02_formatted_data/smallest_demo \
	--output_dir ./data/CTCF_demo/05_motif_discovery
