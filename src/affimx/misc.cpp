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
void print_progress_range()
{
	int i;
	for( i = 0; i < 100; i ++ )
		cout << "_";
	cout << endl;
}

void update_progress_report(
	int *prev_percent,
	double total_count,
	double *progress )
{
	// update the progress report
	(*progress) ++;
	int percent = (*progress) / total_count * 100;
	for( ; (*prev_percent) < percent; (*prev_percent) ++ )
	{
		cout << "^"; cout.flush();
	}
}


///////////////////////////////////////////////////////////////////////////////////////////
bool convert_PWMs_to_PSAMs(
	s_motif *motifs[],
	int num_motifs,
	int convert_to_pfm )
{
	int i, n, x;
	
	// examine all motifs
	for( i = 0; i < num_motifs; i ++ )
		for( x = 0; x < motifs[ i ] ->PWM_width; x ++ ) // examin each position of the motif
		{
			// cout << i << ":" << x << endl;

			// add pseudocounts, and also find the most probable nucleotide of this position
			double max = -1; // the maximum frequency in this position;
			for( n = 0; n < NUM_N_LETTERS; n ++ )
			{
				// the motif should be converted from PWM to PFM
				if( convert_to_pfm )
					motifs[ i ] ->PWM[ n ][ x ] = pow( 10, motifs[ i ] ->PWM[ n ][ x ] );
					
				if( motifs[ i ] ->PWM[ n ][ x ] >= 0 ) // this frequency is valid
				{
					if( !convert_to_pfm )
						motifs[ i ] ->PWM[ n ][ x ] += PSEUDOCOUNT; // add pseudocount
					// NOTE: pseudocount is added only if the input is PFM, not if the input is PWM that is converted to PWM
					
					if( n == 0 || max < motifs[ i ] ->PWM[ n ][ x ] )
					// this is the maximum frequency found so far
						max = motifs[ i ] ->PWM[ n ][ x ];
				}
				else // invalid frequency
				{
					cout << "ERROR: The PWM " << motifs[ i ] ->motif_name
						<< " contains a negative value for nucleotide "
						<< __n_letters[ n ] << " at position " << x+1 << "." << endl;
						
					return false;
				}
			}
				
			if( max <= 0 ) // all nucleotides have zero frequency
			{
				cout << "ERROR: Zero frequency at position " << x+1 << " of the PWM "
					<< motifs[ i ] ->motif_name << "." << endl;
				return false;
			}
			
			// now, normalize the frequencies, so that the maximum frequency of this position is 1.0
			for( n = 0; n < NUM_N_LETTERS; n ++ )
				motifs[ i ] ->PWM[ n ][ x ] /= max;
		}

	return true;
}

///////////////////////////////////////////////////////////////////////////////////////////
void shuffle_motifs(
	s_motif *motifs[],
	int num_motifs )
{
	int i, n, x;
	
	for( i = 0; i < num_motifs; i ++ ) 	// shuffle each motif i
	{
		// shuffle columns
		for( x = 0; x < motifs[ i ] ->PWM_width - 1; x ++ ) // swap each column x with a random column
		{
			int rnd_index = ( rand() % (motifs[ i ] ->PWM_width-x) ) + x;

			// swap column x with column rnd_index
			for( n = 0; n < NUM_N_LETTERS; n ++ )
			{
				double swap = motifs[ i ] ->PWM[ n ][ x ];
				motifs[ i ] ->PWM[ n ][ x ] = motifs[ i ] ->PWM[ n ][ rnd_index ];
				motifs[ i ] ->PWM[ n ][ rnd_index ] = swap;
			}				
		}
		
		// within each column, swap the rows
		for( x = 0; x < motifs[ i ] ->PWM_width; x ++ )
			for( n = 0; n < NUM_N_LETTERS - 1; n ++ )
			{
				int rnd_index = ( rand() % (NUM_N_LETTERS-n) ) + n;

				double swap = motifs[ i ] ->PWM[ n ][ x ];
				motifs[ i ] ->PWM[ n ][ x ] = motifs[ i ] ->PWM[ rnd_index ][ x ];
				motifs[ i ] ->PWM[ rnd_index ][ x ] = swap;
			}
	}
				
}

///////////////////////////////////////////////////////////////////////////////////////////
void shuffle_sequences(
	s_gene *genes[],
	int num_genes )
{
	int i, n, x;
	
	for( i = 0; i < num_genes; i ++ ) 	// shuffle each gene i
		for( x = 0; x < genes[ i ] ->seq_length - 1; x ++ ) // swap each position x with a random position
		{
			int range = MIN( SEQ_SHUFFLE_RANGE, genes[ i ] ->seq_length - x );
			//int range = genes[ i ] ->seq_length - x;
			int rnd_index = ( rand() % range ) + x;

			// swap position x with column rnd_index
			char swap = genes[ i ] ->seq[ x ];
			genes[ i ] ->seq[ x ] = genes[ i ] ->seq[ rnd_index ];
			genes[ i ] ->seq[ rnd_index ] = swap;
		}
}

///////////////////////////////////////////////////////////////////////////////////////////
bool _is_delimiter(
	char ch,
	const char *delimiters,
	int num_delimiters )
// determines whether a character is a delimiter or not
// returns true if ch is a delimiter, false otherwise
//
// ch: the query character
// delimiters: the list of delimiters
// num_delimiters: the size of the delimiter list
{
	// examine the delimiters one by one to find a match
	for( int i = 0; i < num_delimiters; i ++ )
		if( ch == delimiters[ i ] )
			return true; // a matching delimiter is found

	return false;
}


///////////////////////////////////////////////////////////////////////////////////////////
char* extract_phrase(
	const char *line,
	int phrase_index,
	char *phrase,
	const char *delimiters )
// extracts a phrase from a delimitted text
// returns NULL if unsuccessful, 'phrase' otherwise
//
// line: the line from which the phrase will be extracted
// phrase_index: the index of the phrase that will be extracted; 0 means that first phrase, and so on
// phrase: the string in which the extracted phrase will be stored
// delimiters: the list of valid delimiters, separating different phrases within the text
{
	int len = strlen( line );
	int num_delimiters = strlen( delimiters );

	// find the first position after 'phrase_index'th occurrance of a delimiter
	int curr_index = 0;
	int start_pos;
	for( start_pos = 0; start_pos < len && curr_index < phrase_index; start_pos ++ )
		if( _is_delimiter( line[ start_pos ], delimiters, num_delimiters ) )
			curr_index ++;

	// return NULL if the requested phrase is not found
	if( start_pos >= len )
		return NULL;

	// extract the phrase
	int pos;
	for( pos = start_pos; pos < len; pos ++ )
		if( _is_delimiter( line[ pos ], delimiters, num_delimiters ) )
			break;
		else
			phrase[ pos - start_pos ] = line[ pos ];

	phrase[ pos - start_pos ] = 0;

	// return the extracted phrase
	return phrase;
}
