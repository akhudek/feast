#include "dot_plot.h"

void dot_plot_coords( std::ostream &out, ex_sequence_ptr_vector_ptr aln, int start_a, int start_b, int id ) {
	using namespace btl;
	using namespace std;
	ex_sequence::data_type &A = (*aln)[0]->data, &B = (*aln)[1]->data;
	unsigned int next_a = (unsigned int) start_a, next_b = (unsigned int) start_b;

	assert(A.size() == B.size() );

	// initialize window
	for( unsigned int i = 0; i < A.size(); ++i ) {
		if (A[i] != mdna_iupac_gap::GAP && B[i] != mdna_iupac_gap::GAP ) {
			out << next_a << '\t' << next_b << '\t' << id << endl;
			next_a++, next_b++;
		} else 	if (A[i] != mdna_iupac_gap::GAP && B[i] == mdna_iupac_gap::GAP ) next_a++;
		else 	if (A[i] == mdna_iupac_gap::GAP && B[i] != mdna_iupac_gap::GAP ) next_b++;
	}
}
