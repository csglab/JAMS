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

//////////////////////////////////////////////////////////
int _n_index( char nucleotide )
// returns the index of the specified nucleotide
{
	int i;
	for( i = 0; i < NUM_N_LETTERS; i ++ )
		if( nucleotide == __n_letters[ i ] )
			return i;
			
	return -1;
}

//////////////////////////////////////////////////////////
double _get_binding_p(
	char *seq,
	s_motif *motif,
	int direction,
	int *max_dir )
// This function returns the score of a particular sequence stretch relative to the given PWM
{
	double score_f = 1;
	double score_r = 1;
	
	int i;
	for( i = 0; i < motif ->PWM_width; i ++ )
	{
		int pos_f = i;
		int pos_r = motif ->PWM_width - i - 1;
		
		double this_score_f = 0;
		double this_score_r = 0;
		
		int n_index = _n_index( seq[ i ] );

		if( n_index < 0 ) // the nucleotide is unknown
		// thus, the scores will be averaged over this position of the PWM
			for( n_index = 0; n_index < NUM_N_LETTERS; n_index ++ )
			{
				this_score_f += motif ->PWM[ n_index ][ pos_f ] / double(NUM_N_LETTERS);
				this_score_r += motif ->PWM[ n_index ][ pos_r ] / double(NUM_N_LETTERS);
			}
		else // the nucleotide is known
		{
			this_score_f = motif ->PWM[ n_index ][ pos_f ];
			this_score_r = motif ->PWM[ NUM_N_LETTERS - n_index - 1 ][ pos_r ];
		}
		
		// since we're actually dealing with PSAMs, the probabilities are multiplied
		score_f *= this_score_f;
			
		// since we're actually dealing with PSAMs, the probabilities are multiplied
		score_r *= this_score_r;
	}

	if( direction == FORWARD_ONLY ) // only the forward direction is requested
	{
		*max_dir = 1;
		return score_f;
	}
	else if( direction == REVERSE_ONLY ) // only the reverse direction is requested
	{
		*max_dir = -1;
		return score_r;
	}

	// any other value for direction means that both directions should be considered
	if( score_f >= score_r )
		*max_dir = 1;
	else
		*max_dir = -1;

	return score_f + score_r;
}

///////////////////////////////////////////////////////////////////////////////////////////
double scan_sequence(
	s_gene *gene,
	int direction,
	s_motif *motif,
	int *max_pos )
// This function scans the given gene for the instances of the motif
//  and returns the sum of the scores, and the position where maximum score occurs
//
{
	// now calculate the sum of the scores for this motif over the specified bin		
	double score = 0;
	double max_p = -1;
	int end_point = gene ->seq_length - motif ->PWM_width; // the bin cannot end after the last base where a potential motif could be fit
	int i;
	for( i = 0; i <= end_point; i ++ )
	{
		int max_dir;
		double binding_p = _get_binding_p( gene ->seq + i, motif, direction, &max_dir );
		binding_p = 1.0 / ( 1.0 + ( 1.0/binding_p ) );
		score += binding_p;
		
		if( i == 0 || max_p < binding_p ) // this is the max score position so far
		{
			max_p = binding_p;
			
			if( max_dir > 0 )
				*max_pos = i+1;
			else
				*max_pos = -(i + motif ->PWM_width);
		}
	}
		
	return score;
}
