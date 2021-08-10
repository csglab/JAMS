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

#ifndef _H_STRUCTURES // this is to make sure that this file will not be included twice
//*******************************************************************
#define _H_STRUCTURES // this macro indicates that this file is now included

#ifndef _NAMESPACE
using namespace std;
#define _NAMSPACE
#endif


#include <string.h>
#include <stdlib.h>

// ************************************************ macro definitions

// constants
#define MAX_GENES				300000
#define MAX_SEQ_LENGTH			3000000
#define MAX_MOTIFS				100000
#define MAX_PWM_WIDTH			1000
#define MAX_LINE_LENGTH			100000
#define MAX_STRING_LENGTH		1000

#define SHUFFLE_NONE			0
#define SHUFFLE_MOTIFS			1
#define SHUFFLE_SEQS			2

#define SEQ_SHUFFLE_RANGE		MAX_SEQ_LENGTH

#define PSEUDOCOUNT				0.00001

#define FORWARD_ONLY			0
#define REVERSE_ONLY			1
#define BOTH_DIRECTIONS			2

// macros
#define RELEASE(x)		do{ if(x) delete (x); }while(false) // this do-while loop is called only once
#define _CALL(x)		do{ if(!(x)) return false; }while(false)
#define _COPY_STR(x,y)	do{ x=new char[strlen(y)+1]; strcpy((x),(y)); }while(false)
#define MAX(x,y)		(((x)>(y))?(x):(y))
#define MIN(x,y)		(((x)<(y))?(x):(y))

// ********************************************** typedef definitions

typedef unsigned char BYTE;

// ********************************************** local variables
static char __n_letters[] = "ACGT";
#define NUM_N_LETTERS	4


// ******************************************** structure definitions

// this structure will hold the information for the input genes
struct s_gene
{
	s_gene()
	{
		// initialize all the variables, setting them to zero
		
		name = NULL; // no memory allocated yet
		index = -1;
		seq = NULL; // no memory allocated yet
		seq_length = 0; // the initial length of the sequence is zero
	}
	~s_gene()
	{
		// release all the allocated memory
	
		RELEASE( name );
		RELEASE( seq );
	}

	char *name; // the name of this gene, as read from the input FASTA file
	int index;

	char *seq; // the gene, as read from the input FASTA file
	int seq_length; // the length of the nucleotide sequence
};

// this structure will hold the information for the motifs
struct s_motif
{

	s_motif()
	{
		gene_name = NULL;
		motif_name = NULL;
		memset( PWM, 0, sizeof(double*) * NUM_N_LETTERS );
		PWM_width = 0;
	}
	~s_motif()
	{
		RELEASE( gene_name );
		RELEASE( motif_name );
		
		int n;
		for( n = 0; n < NUM_N_LETTERS; n ++ )
			RELEASE( PWM[ n ] );
	}

	char *gene_name;
	char *motif_name;
	double *PWM[ NUM_N_LETTERS ];
	int PWM_width;
};

//*******************************************************************
#endif // this is to make sure that this file will not be included twice
