
#include "annotation_functions.h"


void write_seqa_gff_annotation( std::ostream &out, annotation const &a, pairwise_dna_alignment &align ) {
	using namespace std;
	int col_a = -1, pos_a = -1, first_pos = -1;
	for( annotation_citer i = a.begin(); i != a.end(); ++i ) {
		// find the first position in the first sequence 
		first_pos = -1;
 		for( int c = i->a; c <= i->b; ++c ) {
			++col_a;
			if( align.a->data[col_a] != dna_alignment_alpha::GAP ) {
				++pos_a;
				if( first_pos < 0 ) first_pos = pos_a;
			}
		}
		// now first_pos is the start in seqa and pos_a is the end

		out << align.a->tags["accession"] << "\tCAPE\t";
		if( i->t == 0 ) out << "random\t";
		else out << "model" << i->t << "\t";
		
		out << (first_pos+1) << "\t" << (pos_a+1) << "\t";
		out << ".\t.\t.\t." << endl;
	}
}
