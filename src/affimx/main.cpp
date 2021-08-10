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


// main.cpp : Defines the entry point for the console application.
//

#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <math.h>
#include <time.h>

#include "declarations.h"

extern const char *__PWM_file;
extern const char *__FASTA_file;
extern const int *__direction;
extern const int *__shuffle;
extern const char *__output_file;

unsigned int _strtoint( const char *string )
{
	unsigned int rval = 0;
	int len = strlen( string );
	int i;
	for( i = 0; i < len; i ++ )
		rval += string[ i ];
		
	return rval;
}
		
void welcome_message()
{
	cout << endl
		<< "***************************** AffiMx version 9.0 ******************************" << endl
		<< "*                                                                             *" << endl
		<< "* Copyright 2016 Hamed S. Najafabadi                                          *" << endl
		<< "*                                                                             *" << endl
		<< "*******************************************************************************" << endl
		<< endl;
}


int main( int argc, char* argv[] )
{	
	welcome_message();

	if( argc <= 1 )
	// no argument is profided
	{
		print_commandline_format();
		return 0;
	}
	
	if( !read_arguments( argc, argv ) )
		return 1;
		
	print_arguments();

	//******************* open output
	
	ofstream ofs_aff;
	if( !open_output( ofs_aff, __output_file, ".affinity.txt" ) )
		return 1;

	ofstream ofs_ene;
	if( !open_output( ofs_ene, __output_file, ".energy.txt" ) )
		return 1;

	ofstream ofs_pos;
	if( !open_output( ofs_pos, __output_file, ".position.txt" ) )
		return 1;

	//******************* open the FASTA file

	s_gene *genes[ MAX_GENES ];
	int num_genes = 0;

	cout << "Opening the FASTA file..." << endl;
	if( !open_FASTA( __FASTA_file, genes, &num_genes ) )
		return 1;
	
	//******************* open the PWM file

	s_motif *motifs[ MAX_MOTIFS ];
	int num_motifs = 0;

	int is_PWM = 0;
	cout << "Opening the PWM file..." << endl;
	if( !open_motifs( __PWM_file, motifs, &num_motifs, &is_PWM ) )
		return 1;
		
	cout << "Converting the matrices to PSAMs..." << endl;
	if( !convert_PWMs_to_PSAMs( motifs, num_motifs, is_PWM ) )
		return 1;
		
	// set the srand value to a number that depends on time
	srand( time( NULL ) + _strtoint( __output_file ) );
	if( *__shuffle == SHUFFLE_MOTIFS )
	{
		cout << "Randomizing the motifs..." << endl;
		shuffle_motifs( motifs, num_motifs );
	}
	else if( *__shuffle == SHUFFLE_SEQS )
	{
		cout << "Randomizing the sequences..." << endl;
		shuffle_sequences( genes, num_genes );
	}

	//******************* write the output
	cout << "Writing the outputs..." << endl;
	write_affinity_matrix( ofs_aff, ofs_ene, ofs_pos,
		genes, num_genes, motifs, num_motifs,
		*__direction );
	
	cout << endl << "Job finished successfully." << endl;

	return 0;
}
