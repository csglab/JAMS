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


#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <math.h>
#include <time.h>

#include "declarations.h"

///////////////////////////////////////////////////////////////////////////////////////////
void write_affinity_matrix(
	ofstream &ofs_aff,
	ofstream &ofs_ene,
	ofstream &ofs_pos,
	s_gene *genes[],
	int num_genes,
	s_motif *motifs[],
	int num_motifs,
	int direction )
// writes the affinity matrix
//
// ofs_aff: the output file stream for the affinities
// ofs_ene: the output file stream for the energies (log of affinities)
// genes: the gene list
// num_genes: the gene list size
// motifs: the motif list
// num_motifs: the motif list size
// direction: see structures.h
{
	// variables for progress report
	int prev_percent = 0;
	double total_count = num_genes;
	double progress = 0;
	print_progress_range();
		
	//////////////////////////////////////////////// print the header
	ofs_aff << "Name";
	ofs_ene << "Name";
	ofs_pos << "Name";

	int i, j, k;
	for( i = 0; i < num_motifs; i ++ )
	{
		ofs_aff << char(9) << motifs[ i ] ->motif_name;
		ofs_ene << char(9) << motifs[ i ] ->motif_name;
		ofs_pos << char(9) << motifs[ i ] ->motif_name;
	}
	ofs_aff << endl;
	ofs_ene << endl;
	ofs_pos << endl;
	
	for( i = 0; i < num_genes; i ++ )
	{
		ofs_aff << genes[ i ] ->name;
		ofs_ene << genes[ i ] ->name;
		ofs_pos << genes[ i ] ->name;
		
		int j;
		for( j = 0; j < num_motifs; j ++ )
		{
			int max_pos = 0;
			double sum_score = scan_sequence( genes[ i ], direction, motifs[ j ], &max_pos );

			ofs_aff << char(9) << sum_score;
			ofs_ene << char(9) << log10( sum_score );
			ofs_pos << char(9) << max_pos;
		}

		ofs_aff << endl;
		ofs_ene << endl;
		ofs_pos << endl;

		// update the progress report
		update_progress_report( &prev_percent, total_count, &progress );
	}
	
	cout << endl;
}
