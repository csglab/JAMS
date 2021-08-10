// Copyright 2012 Hamed Shateri Najafabadi

/********************************************************************

This file is part of AffiMx.

AffiMx is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

AffiMx is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with AffiMx.  If not, see <http://www.gnu.org/licenses/>.

********************************************************************/

#ifndef _H_DECLARATIONS // this is to make sure that this file will not be included twice
//*******************************************************************
#define _H_DECLARATIONS // this macro indicates that this file is now included


#include "structures.h" // this header file is required for definition of structures


bool read_arguments( int argc, char* argv[] );
void print_arguments();
void print_commandline_format();

void print_progress_range();
void update_progress_report( int *prev_percent, double total_count, double *progress );
bool convert_PWMs_to_PSAMs( s_motif *motifs[], int num_motifs, int convert_to_pfm );
void shuffle_motifs( s_motif *motifs[], int num_motifs );
void shuffle_sequences( s_gene *genes[], int num_genes );
char* extract_phrase( const char *line, int phrase_index, char *phrase, const char *delimiters );

bool open_output( ofstream &ofs, const char *path, const char *extension );
bool open_FASTA( const char *FASTA_file, s_gene *genes[], int *num_genes );
bool open_motifs( const char *motif_file, s_motif *motifs[], int *num_motifs, int *is_PWM );

bool read_FASTA( ifstream &ifs, s_gene *genes[], int *num_genes, char this_EOL );
bool read_motifs( ifstream &ifs, s_motif *motifs[], int *num_motifs, int *is_PWM, char this_EOL );

double scan_sequence( s_gene *gene, int direction, s_motif *motif, int *max_pos );

void write_affinity_matrix( ofstream &ofs_aff, ofstream &ofs_ene, ofstream &ofs_pos, s_gene *genes[], int num_genes, s_motif *motifs[], int num_motifs, int direction );

//*******************************************************************
#endif // this is to make sure that this file will not be included twice
