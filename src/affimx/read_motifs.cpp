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

static int _line_number; // this local variable will hold the number of lines that are read so far
	
static char _delimiters[] = { char(9), 0 };

/////////////////////////////////////////////////////////////////////////////////
bool _add_motif(
	char gene_name[],
	char motif_name[],
	double *PWM[],
	int PWM_width,
	s_motif *motifs[],
	int *num_motifs )
// adds an entry to the gene list
//
// gene_name: the name of the new gene
{
	if(	!PWM_width ) // there is no associated PWM, and thus no need to do anything
		return true; // return with no error
		
	if( !motif_name[ 0 ] ) // there is an associated PWM, but no motif name
	{
		cout << "ERROR: No motif name provided (error at line " << _line_number << ")"
			<< endl;
		return false;
	}
	
	if( *num_motifs >= MAX_MOTIFS ) // too many motifs
	{
		cout << "ERROR: Maximum number of motifs (" << MAX_MOTIFS << ") reached "
			<< "(error at line " << _line_number << ")" << endl;
		return false;
	}
	
	// create the motif
	motifs[ *num_motifs ] = new s_motif;
	
	// store the gene and motif names
	_COPY_STR( motifs[ *num_motifs ] ->gene_name, gene_name );
	_COPY_STR( motifs[ *num_motifs ] ->motif_name, motif_name );
	
	// create and store the PWM
	int n;
	for( n = 0; n < NUM_N_LETTERS; n ++ )
	{
		motifs[ *num_motifs ] ->PWM[ n ] = new double[ PWM_width ]; // allocate memory
		memcpy( motifs[ *num_motifs ] ->PWM[ n ], PWM[ n ], sizeof(double) * PWM_width ); // copy
	}
	
	motifs[ *num_motifs ] ->PWM_width = PWM_width;

	// update the number of motifs	
	(*num_motifs) ++;
	
	return true;
}


/////////////////////////////////////////////////////////////////////////////////
bool read_motifs(
	ifstream &ifs,
	s_motif *motifs[],
	int *num_motifs,
	int *is_PWM,
	char this_EOL )
// reads motifs from &ifs
// returns true if successful, false otherwise
//
// ifs: input file stream
// motifs: the list of motifs whose PWMs are to be read
// num_motifs: will contain the number of motifs that are in the list
// this_EOL: the end of line character. this character is determined separately in another function for this file stream
{
	*is_PWM = 0;

	// A line cannot be longer than MAX_LINE_LENGTH characters
	char string[ MAX_LINE_LENGTH + 1 ];
	
	char gene_name[ MAX_LINE_LENGTH + 1 ] = ""; // not used in this version
	char motif_name[ MAX_LINE_LENGTH + 1 ] = "";
	double *PWM[ NUM_N_LETTERS ];
	int PWM_width = 0;
	// initialize the PWM memory
	int n;
	for( n = 0; n < NUM_N_LETTERS; n ++ )
	{
		PWM[ n ] = new double[ MAX_PWM_WIDTH ];
		memset( PWM[ n ], 0, sizeof(double) * MAX_PWM_WIDTH );
	}


	///////////// read the file line by line

	_line_number = 0; // reset the line number
	while( true )
	{
		// read the line
		ifs.getline( string, MAX_LINE_LENGTH, this_EOL );

		if( !ifs || ifs.eof() ) // the end of the file has reached
		{
			// add the last read motif to the list
			_CALL(
				_add_motif( gene_name, motif_name, PWM, PWM_width, motifs, num_motifs )
				);

			break;
		}
			
		_line_number ++; // update the line number that was just read

		char first_phrase[ MAX_STRING_LENGTH + 1 ];
		if( !extract_phrase( string, 0, first_phrase, _delimiters ) )
			continue; // empty line
		
		if( strcmp( first_phrase, "Motif" ) == 0 ) // this line indicates the gene name
		{
			// add the last read motif to the list
			_CALL(
				_add_motif( gene_name, motif_name, PWM, PWM_width, motifs, num_motifs )
				);

			// this is a new motif
			motif_name[ 0 ]= 0;
			PWM_width = 0;

			if( !extract_phrase( string, 1, motif_name, _delimiters ) )
			{
				cout << "ERROR: Motif name expected at line " << _line_number << endl;
				return false;
			}
			
		}
		else if( strchr( "0123456789", first_phrase[0] ) )
		// this is a matrix line
		{
			int pos = atoi( first_phrase ) - 1;
			PWM_width = pos + 1; // the length of the PWM is the index of the last position plus one
			
			for( n = 0; n < NUM_N_LETTERS; n ++ )
			{
				char str_value[ MAX_STRING_LENGTH + 1 ];
				if( !extract_phrase( string, n + 1, str_value, _delimiters ) )
				{
					cout << "ERROR: Unexpected end of line at line " << _line_number << endl;
					return false;
				}
				
				PWM[ n ][ pos ] = atof( str_value );
				
				if( PWM[ n ][ pos ] < 0 )
					*is_PWM = 1;
			}
		}
	}
	
	return true;
}
