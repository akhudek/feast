/** Generic Pair-HMM implementation.
 *
 *  This class is poorly organized at the moment. Constrained Baum-Welch
 *  training is implemented, but not Viterbi. Various flavours of seed
 *  extensions are also implemented here.
 */

#ifndef PAIR_HMM_H__
#define PAIR_HMM_H__

#include "ereal.h"

#include <vector>
#include <set>
#include <blitz/array.h>
#include <boost/tuple/tuple.hpp>
#include <boost/shared_ptr.hpp>
#include <boost/thread.hpp>
#include "extension_store.h"
#include "closed_interval.h"
#include "concurrent_queue.h"
#include "point.h"

template<typename ALPHA, typename DT = ereal >
class pair_hmm;

template<typename ALPHA, typename DT, typename D1, typename D2, typename ITERATOR >
class pair_hmm_forward_worker {
protected:
	typedef blitz::Array<DT,2> border;
	typedef blitz::Array<DT,3> data;

	pair_hmm<ALPHA,DT> const &master;

	ITERATOR start;

	ITERATOR end;

	int interface_size;

	std::vector<border> &forward_borders;

public:
	pair_hmm_forward_worker( pair_hmm<ALPHA,DT> const &m, ITERATOR st, ITERATOR en, int isize, std::vector<border> &fb )
	: master(m), start(st), end(en), interface_size(isize), forward_borders(fb) {}

	void operator()() {
		using namespace std;
		using namespace blitz;
		Range all = Range::all();

		D1 firstRa = *(start->first);
		border first_border(firstRa.data.size()+1,master.M), nullbd;
		master.init_forward( border(), first_border, firstRa.data);
		forward_borders.push_back( first_border(Range(0,firstRa.data.size()-1),all).copy() );

		data f;

		ITERATOR next = start, it = start; ++next;

		while( 1 ) {
			// first adjust region to align based on interface
			D1 rega = *(it->first);	D2 regb = *(it->second);

//			assert( rega.data.size() >= interface_size );

			// if not at start, add upper interface. Previous region included "interface" bases from a.
			if( it != start ) rega.data.range.a -= interface_size;

			// if not at end, add lower interface
			if( next != end ) rega.data.range.b += interface_size;

			// allocate dp matrices and fill
			int n = (int) rega.data.size();
			f.resize(n+1,2,master.M); f = 0.0;

			master.init_forward(forward_borders.back(),f(all,0,all),rega.data);
			//if( status ) std::cerr << "Doing forward computation " << rega.data.range << " vs " << regb.data.range << std::endl;

			// compute breakpoints and final forward score
			int col = master.compute_breakpoints( 0, data(), f, rega.data, regb.data );

			// Copy border data back to idata. We actually store interface_size + 1 in this case, since we
			// initialize the border element too.
			if( next == end ) forward_borders.push_back( f(Range(n-1,n),col,all).copy() );
			else forward_borders.push_back( f(Range(n-interface_size*2,n),col,all).copy() );

			if( next == end ) break;

			++next;
			++it;

		}
	}
};


template<typename ALPHA, typename DT, typename D1, typename D2, typename ITERATOR >
class pair_hmm_backward_worker {
protected:
	typedef blitz::Array<DT,2> border;
	typedef blitz::Array<DT,3> data;

	pair_hmm<ALPHA,DT> const &master;

	ITERATOR start;

	ITERATOR end;

	int interface_size;

	std::queue<border> &backward_borders;

public:
	pair_hmm_backward_worker( pair_hmm<ALPHA,DT> const &m, ITERATOR st, ITERATOR en, int isize, std::queue<border> &bb )
	: master(m), start(st), end(en), interface_size(isize), backward_borders(bb) {}

	void operator()() {
		using namespace blitz;
		Range all = Range::all();
		ITERATOR it = end, next = end; --it;

		D1 lastRa = *(it->first);
		if( it != start ) lastRa.data.range.a -= interface_size;
		border backward_border_init(lastRa.data.size()+1,master.M);
		master.init_backward( border(), backward_border_init, lastRa.data );
		backward_borders.push( backward_border_init(Range(std::max(1,(int)lastRa.data.size()-interface_size*2),lastRa.data.size()),all).copy() );

		data b;

		while( it != start ) {
			// setup region
			D1 rega = *(it->first);	D2 regb = *(it->second);

			assert( rega.data.size() >= (unsigned int)interface_size );

			// we are never at start
			rega.data.range.a -= interface_size;
			regb.data.range.a -= 1;

			// if not at end, add interface at end
			if( next != end ) rega.data.range.b += interface_size;

			// allocate dp matrices and fill
			int n = (int) rega.data.size();
			b.resize(n+1,2,master.M); b = 0.0;

			master.init_backward(backward_borders.back(),b(all,0,all),rega.data);
			int col = master.compute_backward( b, rega.data, regb.data );

			// copy border data back to idata
			backward_borders.push( b(Range(1,interface_size*2),col,all).copy() );

			--it;
			--next;
		}
	}
};


template<typename ALPHA, typename DT, typename D1, typename D2, typename ITERATOR >
class pair_hmm_bw_worker {
	typedef blitz::Array<DT,2> border;
	typedef blitz::Array<DT,3> data;

public:

	struct worker_data {
		border fw_border;
		border bw_border;
		D1 A;
		D2 B;
		worker_data( border &fw, border &bw, D1 a, D2 b ) : A(a), B(b) {
			fw_border.resize( fw.shape() );
			fw_border = fw;
			bw_border.resize( bw.shape() );
			bw_border = bw;
		}
		worker_data() {}

		worker_data &operator=( worker_data const &o ) {
			fw_border.resize( o.fw_border.shape() );
			fw_border = o.fw_border;
			bw_border.resize( o.bw_border.shape() );
			bw_border = o.bw_border;
			A = o.A;
			B = o.B;
			return *this;
		}
	};

	typedef concurrent_queue<worker_data> work_queue;
protected:

	pair_hmm<ALPHA,DT> const &master;

	work_queue &wqueue;

	blitz::Array<DT,2> Et;
	blitz::Array<DT,3> Ee;

public:
	pair_hmm_bw_worker( pair_hmm<ALPHA,DT> const &m, work_queue &wq, blitz::Array<DT,2> T, blitz::Array<DT,3> E )
	: master(m), wqueue(wq), Et(T), Ee(E) {}

	void operator()() {
		worker_data wd;
		while( wqueue.get_data(wd) ) {
			//Et = 0.0; Ee = 0.0;
			master.train_bw_checkpoint( wd.fw_border, wd.bw_border, wd.A.data, wd.B.data, Et, Ee );
			//master.train_bw( wd.fw_border, wd.bw_border, wd.A.data, wd.B.data, Et, Ee );

			//std::cerr << Ee(4,blitz::Range::all(),blitz::Range::all()) << std::endl;
		}
	}
};

/** Pair-HMM class.
 * ALPHA is the alphabet. DT is the storage type and defaults to
 * ereal.
 */
template<typename ALPHA, typename DT  >
class pair_hmm {
public:
	/** Symbol indicating no output.
	 * This is always the last symbol in the alphabet.
	 */
	static int const NOTHING = ALPHA::SIZE;

	/** Number of states. */
	int M;

	/** Transition pr: T(i,j) - from state i to j. */
	blitz::Array<DT,2> T;

	/** Emmision pr: E(m,a,b) - emit a in seq 1 and b in seq 2 from state m. */
	blitz::Array<DT,3> E;

	/** Start state. */
	int start_state;

	/** End state. */
	int end_state;

	// order in which to evaluate states
	std::vector<int> eval_order, eval_orderB;

	// dependancy lists
	blitz::Array< std::vector<int>, 1 > depS, depBS;
	blitz::Array< std::vector<int>, 2 > depL, depR, depU, depD, depAB;
	blitz::Array< std::vector<int>, 3 > depUL, depDR;

	// operation list
	enum op_direction { MM = 0, MZ, ZM, ZZ, ZERO };
	typedef boost::tuple< int,unsigned char,int,ereal, bool > operation;
	typedef std::vector<operation> op_vector;
	typedef std::vector< std::vector<op_vector> > op_table;
	op_table operations;

	blitz::Array <std::vector<operation>, 2> op_list;

	// Data for forward algorithm.
	typedef blitz::Array<DT,1> cell;		// d(i,cell)
	typedef blitz::Array<DT,2> border;		// d(i,state)
	typedef blitz::Array<DT,3> data;		// d(i,j,state)

	// Stores a path through the pair hmm as a compressed tree.
	typedef boost::tuple<int,int,int> path_column;
	typedef std::vector<path_column> simple_path;
	struct compressed_path {
		simple_path data;
		boost::shared_ptr<compressed_path> chained_path;
		int chained_pos;
	};

	class path {
	protected:
		boost::shared_ptr<compressed_path> cp;
		int end_pos;

	public:
		path() : cp( new compressed_path ), end_pos(0) {}

		/** Set path to be path o plus column c. */
		void set( path &o, path_column c ) {
			// Path o holds the end of it, so we extend c.
			if( o.cp->data.size() == o.end_pos ) {
				cp = o.cp;
			// Path o has been extended by some other path, so we need to branch.
			} else {
				boost::shared_ptr<compressed_path> mycp( new compressed_path );
				cp = mycp;
				cp->chained_path = o.cp;
				cp->chained_pos  = o.end_pos;
			}
			cp->data.push_back(c);
			end_pos = cp->data.size();
		}

		void copy_to_simple_path( simple_path &sp ) {
			using namespace std;
			vector< pair< boost::shared_ptr<compressed_path>,int > > stack;
			stack.push_back( make_pair(cp,end_pos) );
			while( stack.back().first->chained_path ) stack.push_back( make_pair(stack.back().first->chained_path,stack.back().first->chained_pos) );

			while( ! stack.empty() ) {
				copy( stack.back().first->data.begin(), stack.back().first->data.begin() + stack.back().second, back_inserter(sp) );
				stack.pop_back();
			}
		}
	};

	template< typename D1 >
	void init_forward( border bd, border f, D1 A ) const {
		using namespace blitz;
		using namespace std;
		int n = (int) A.size();
		Range all = Range::all();

		assert( f.extent(firstDim) >= n+1 );
		assert( f.extent(firstDim) >= bd.extent(firstDim) );

		f(all,all) = 0.0;

		// initialize with border data or just forward border data, not both
		if( bd.extent(firstDim) > 0 ) f(Range(0,bd.extent(firstDim)-1),all) = bd;
		int i_start = bd.extent(firstDim);
		if( i_start == 0 ) {
			f(0,start_state) = 1.0;
			for( vector<int>::const_iterator k = eval_order.begin(); k != eval_order.end(); ++k )
				for( vector<int>::const_iterator s = depS(*k).begin(); s != depS(*k).end(); ++s )
					f(0,*k) += f(0,*s)*T(*s,*k)*E(*k,NOTHING,NOTHING);
			i_start = 1;
		}

		for( int i = i_start; i <= n; ++i ) {
			// forward probability
			for( vector<int>::const_iterator k = eval_order.begin(); k != eval_order.end(); ++k ) {
				for( vector<int>::const_iterator u = depU(*k,(int)A[i-1]).begin(); u != depU(*k,(int)A[i-1]).end(); ++u )
					f(i,*k) += f(i-1,*u)*T(*u,*k)*E(*k,(int)A[i-1],NOTHING);
				for( vector<int>::const_iterator s = depS(*k).begin(); s != depS(*k).end(); ++s )
					f(i,*k) += f(i,*s)*T(*s,*k)*E(*k,NOTHING,NOTHING);
			}
		}
	}

	template< typename D1 >
	void init_backward( border bd, border b, D1 A ) const {
		using namespace blitz;
		using namespace std;
		int n = (int) A.size();
		Range all = Range::all();

		assert( b.extent(firstDim) >= n+1 );
		assert( b.extent(firstDim) >= bd.extent(firstDim) );

		b = 0.0;

		// initialize with border data or just backward border data, not both
		if( bd.extent(firstDim) > 0 ) b(Range(n-bd.extent(firstDim)+1,n),all) = bd;

		int i_start = n-bd.extent(firstDim);
		if( i_start == n ) {
			for( int k = 0; k < M; ++k ) b(n,k) = T(k,end_state);
			for( vector<int>::const_iterator k = eval_orderB.begin(); k != eval_orderB.end(); ++k )
				for( vector<int>::const_iterator s = depBS(*k).begin(); s != depBS(*k).end(); ++s )
					b(n,*k) += T(*k,*s)*b(n,*s)*E(*s,NOTHING,NOTHING);
			i_start = n-1;
		}

		for( int i = i_start; i >= 1; --i ) {
			for( vector<int>::const_iterator k = eval_orderB.begin(); k != eval_orderB.end(); ++k ) {
				for( vector<int>::const_iterator u = depD(*k,(int)A[i]).begin(); u != depD(*k,(int)A[i]).end(); ++u )
					b(i,*k) += T(*k,*u)*E(*u,(int)A[i],NOTHING)*b(i+1,*u);
				for( vector<int>::const_iterator s = depBS(*k).begin(); s != depBS(*k).end(); ++s )
					b(i,*k) += T(*k,*s)*E(*s,NOTHING,NOTHING)*b(i,*s);
			}
		}
	}

	// count the border values
	template< typename D1 >
	void tally_backborder( border bd, border f, border b, D1 A, blitz::Array<DT,2> &Et, blitz::Array<DT,3> &Ee ) const {
		using namespace blitz;
		using namespace std;
		int n = (int) A.size();
		Range all = Range::all();

		assert( b.extent(firstDim) >= n+1 );
		assert( f.extent(firstDim) >= n+1 );

		if( bd.extent(firstDim) == 0 ) {
			for( int k = 0; k < M; ++k ) b(n,k) = T(k,end_state);
			for( vector<int>::const_iterator k = eval_orderB.begin(); k != eval_orderB.end(); ++k )
				for( vector<int>::const_iterator s = depBS(*k).begin(); s != depBS(*k).end(); ++s )
					Et(*k,*s) += f(n,*k)*T(*k,*s)*b(n,*s)*E(*s,NOTHING,NOTHING);

			for( vector<int>::const_iterator s = depAB(NOTHING,NOTHING).begin(); s != depAB(NOTHING,NOTHING).end(); ++s )
				Ee(*s,NOTHING,NOTHING) += f(n,*s)*b(n,*s);
		}

		for( int i = n-bd.extent(firstDim)-1; i > 0; --i ) {
			// forward probability
			for( vector<int>::const_iterator k = eval_orderB.begin(); k != eval_orderB.end(); ++k ) {
				for( vector<int>::const_iterator u = depD(*k,(int)A[i]).begin(); u != depD(*k,(int)A[i]).end(); ++u )
					Et(*k,*u) += f(i,*k)*b(i+1,*u)*T(*k,*u)*E(*u,(int)A[i],NOTHING);
				for( vector<int>::const_iterator s = depBS(*k).begin(); s != depBS(*k).end(); ++s )
					Et(*k,*s) += f(n,*k)*b(i,*s)*T(*k,*s)*E(*s,NOTHING,NOTHING);
			}

			// tally emmissions
			for( vector<int>::const_iterator s = depAB((int)A[i-1],NOTHING).begin(); s != depAB((int)A[i-1],NOTHING).end(); ++s )
				Ee(*s,(int)A[i-1],NOTHING) += f(i,*s)*b(i,*s);
			for( vector<int>::const_iterator s = depAB(NOTHING,NOTHING).begin(); s != depAB(NOTHING,NOTHING).end(); ++s )
				Ee(*s,NOTHING,NOTHING) += f(i,*s)*b(i,*s);
		}

	}

	template< typename D1, typename D2>
	int compute_backward( data b, D1 A, D2 B ) const {
		using namespace std;
		using namespace blitz;
		int n = (int) A.size(), m = (int) B.size();
		Range all = Range::all();

		int prev = 0, cur = 1;
		for( int j = m-1; j >= 1; --j ) {
			// init column
			b(all,cur,all) = 0.0;

			// for i == n we can only do L and S transitions
			for( vector<int>::const_iterator k = eval_orderB.begin(); k != eval_orderB.end(); ++k ) {
				for( vector<int>::const_iterator l = depR(*k,(int)B[j]).begin(); l != depR(*k,(int)B[j]).end(); ++l )
					b(n,cur,*k) += T(*k,*l)*E(*l,NOTHING,(int)B[j])*b(n,prev,*l);
				for( vector<int>::const_iterator s = depBS(*k).begin(); s != depBS(*k).end(); ++s )
					b(n,cur,*k) += T(*k,*s)*E(*s,NOTHING,NOTHING)*b(n,cur,*s);
			}

			for( int i = n-1; i >= 1; --i ) {
				for( vector<int>::const_iterator k = eval_orderB.begin(); k != eval_orderB.end(); ++k ) {
					for( vector<int>::const_iterator r = depR(*k,(int)B[j]).begin(); r != depR(*k,(int)B[j]).end(); ++r )
						b(i,cur,*k) += T(*k,*r)*E(*r,NOTHING,(int)B[j])*b(i,prev,*r);
					for( vector<int>::const_iterator d = depDR(*k,(int)A[i],(int)B[j]).begin(); d != depDR(*k,(int)A[i],(int)B[j]).end(); ++d )
						b(i,cur,*k) += T(*k,*d)*E(*d,(int)A[i],(int)B[j])*b(i+1,prev,*d);
					for( vector<int>::const_iterator u = depD(*k,(int)A[i]).begin(); u != depD(*k,(int)A[i]).end(); ++u )
						b(i,cur,*k) += T(*k,*u)*E(*u,(int)A[i],NOTHING)*b(i+1,cur,*u);
					for( vector<int>::const_iterator s = depBS(*k).begin(); s != depBS(*k).end(); ++s )
						b(i,cur,*k) += T(*k,*s)*E(*s,NOTHING,NOTHING)*b(i,cur,*s);
				}
			}

			// swap columns
			int tmp = cur; cur = prev; prev = tmp;
		}

		return prev;
	}


	template< typename D1, typename D2>
	int compute_breakpoints( int bpspace, data bp, data f, D1 A, D2 B ) const {
		using namespace std;
		using namespace blitz;
		int n = (int) A.size(), m = (int) B.size();
		Range all = Range::all();

		if( bpspace > 0 )	bp(all,0,all) = f(all,0,all);

		int prev = 0, cur = 1;
		for( int j = 1; j <= m; ++j ) {
			// set everything back to zero
			f(all,cur,all) = 0.0;

			// for i == 0 we can only do L and S transitions
			for( vector<int>::const_iterator k = eval_order.begin(); k != eval_order.end(); ++k ) {
				for( vector<int>::const_iterator l = depL(*k,(int)B[j-1]).begin(); l != depL(*k,(int)B[j-1]).end(); ++l )
					f(0,cur,*k) += f(0,prev,*l)*T(*l,*k)*E(*k,NOTHING,(int)B[j-1]);
				for( vector<int>::const_iterator s = depS(*k).begin(); s != depS(*k).end(); ++s )
					f(0,cur,*k) += f(0,cur,*s)*T(*s,*k)*E(*k,NOTHING,NOTHING);
			}

			// now we cover the remaining blocks
			for( int i = 1; i <= n; ++i ) {
				for( vector<int>::const_iterator k = eval_order.begin(); k != eval_order.end(); ++k ) {
					for( vector<int>::const_iterator l = depL(*k,(int)B[j-1]).begin(); l != depL(*k,(int)B[j-1]).end(); ++l )
						f(i,cur,*k) += f(i,prev,*l)*T(*l,*k)*E(*k,NOTHING,(int)B[j-1]);
					for( vector<int>::const_iterator d = depUL(*k,(int)A[i-1],(int)B[j-1]).begin(); d != depUL(*k,(int)A[i-1],(int)B[j-1]).end(); ++d )
						f(i,cur,*k) += f(i-1,prev,*d)*T(*d,*k)*E(*k,(int)A[i-1],(int)B[j-1]);
					for( vector<int>::const_iterator u = depU(*k,(int)A[i-1]).begin(); u != depU(*k,(int)A[i-1]).end(); ++u )
						f(i,cur,*k) += f(i-1,cur,*u)*T(*u,*k)*E(*k,(int)A[i-1],NOTHING);
					for( vector<int>::const_iterator s = depS(*k).begin(); s != depS(*k).end(); ++s )
						f(i,cur,*k) += f(i,cur,*s)*T(*s,*k)*E(*k,NOTHING,NOTHING);
				}
			}

			// save breakpoint
			if( bpspace > 0  && j % bpspace == 0 ) bp(all,j/bpspace,all) = f(all,cur,all);

			// swap columns
			int tmp = cur; cur = prev; prev = tmp;
		}

		return prev;
	}

	// We assume f and b are each of size (n+1,m+1) and that f(0,all), b(n,all)
	// are appropriately initialized.The rest of f should be zero.
	template< typename D1, typename D2 >
	void compute_bw( data f, data b, D1 A, D2 B, blitz::Array<DT,2> &Et, blitz::Array<DT,3> &Ee, bool toend = false ) const {
		using namespace std;
		using namespace blitz;
		int n = (int) A.size(), m = B.size();
		Range all = Range::all();
		assert( f.extent(firstDim) >= n+1 );
		assert( b.extent(firstDim) >= n+1 );
		assert( f.extent(secondDim) >= m+1 );
		assert( b.extent(secondDim) >= m+1 );

		// fill in forward values
		for( int j = 1; j <= m; ++j ) {
			// for i == 0 we can only do L and S transitions
			for( vector<int>::const_iterator k = eval_order.begin(); k != eval_order.end(); ++k ) {
				for( vector<int>::const_iterator l = depL(*k,(int)B[j-1]).begin(); l != depL(*k,(int)B[j-1]).end(); ++l )
					f(0,j,*k) += f(0,j-1,*l)*T(*l,*k)*E(*k,NOTHING,(int)B[j-1]);
				for( vector<int>::const_iterator s = depS(*k).begin(); s != depS(*k).end(); ++s )
					f(0,j,*k) += f(0,j,*s)*T(*s,*k)*E(*k,NOTHING,NOTHING);
			}

			for( int i = 1; i <= n; ++i ) {
				for( vector<int>::const_iterator k = eval_order.begin(); k != eval_order.end(); ++k ) {
					for( vector<int>::const_iterator l = depL(*k,(int)B[j-1]).begin(); l != depL(*k,(int)B[j-1]).end(); ++l )
						f(i,j,*k) += f(i,j-1,*l)*T(*l,*k)*E(*k,NOTHING,(int)B[j-1]);
					for( vector<int>::const_iterator d = depUL(*k,(int)A[i-1],(int)B[j-1]).begin(); d != depUL(*k,(int)A[i-1],(int)B[j-1]).end(); ++d )
						f(i,j,*k) += f(i-1,j-1,*d)*T(*d,*k)*E(*k,(int)A[i-1],(int)B[j-1]);
					for( vector<int>::const_iterator u = depU(*k,(int)A[i-1]).begin(); u != depU(*k,(int)A[i-1]).end(); ++u )
						f(i,j,*k) += f(i-1,j,*u)*T(*u,*k)*E(*k,(int)A[i-1],NOTHING);
					for( vector<int>::const_iterator s = depS(*k).begin(); s != depS(*k).end(); ++s )
						f(i,j,*k) += f(i,j,*s)*T(*s,*k)*E(*k,NOTHING,NOTHING);
				}
			}
		}

		// fill backwards values and gather estimates
		int j_end = toend ? 0 : 1;
		for( int j = m-1; j >= j_end; --j ) {
			// for i == n we can only do L and S transitions
			for( vector<int>::const_iterator k = eval_orderB.begin(); k != eval_orderB.end(); ++k ) {
				for( vector<int>::const_iterator l = depR(*k,(int)B[j]).begin(); l != depR(*k,(int)B[j]).end(); ++l ) {
					DT p = T(*k,*l)*E(*l,NOTHING,(int)B[j])*b(n,j+1,*l);
					b(n,j,*k) += p;
					Et(*k,*l) += f(n,j,*k)*p;
				}
				for( vector<int>::const_iterator s = depBS(*k).begin(); s != depBS(*k).end(); ++s ) {
					DT p = T(*k,*s)*E(*s,NOTHING,NOTHING)*b(n,j,*s);
					b(n,j,*k) += p;
					Et(*k,*s) += f(n,j,*k)*p;
				}
			}

			for( vector<int>::const_iterator s = depAB(NOTHING,(int)B[j-1]).begin(); s != depAB(NOTHING,(int)B[j-1]).end(); ++s )
				Ee(*s,NOTHING,(int)B[j-1]) += f(n,j,*s)*b(n,j,*s);
			for( vector<int>::const_iterator s = depAB(NOTHING,NOTHING).begin(); s != depAB(NOTHING,NOTHING).end(); ++s )
				Ee(*s,NOTHING,NOTHING) += f(n,j,*s)*b(n,j,*s);

			for( int i = n-1; i >= 1; --i ) {
				for( vector<int>::const_iterator k = eval_orderB.begin(); k != eval_orderB.end(); ++k ) {
					for( vector<int>::const_iterator r = depR(*k,(int)B[j]).begin(); r != depR(*k,(int)B[j]).end(); ++r ) {
						DT p = T(*k,*r)*E(*r,NOTHING,(int)B[j])*b(i,j+1,*r);
						b(i,j,*k) += p;
						Et(*k,*r) += f(i,j,*k)*p;
					}
					for( vector<int>::const_iterator d = depDR(*k,(int)A[i],(int)B[j]).begin(); d != depDR(*k,(int)A[i],(int)B[j]).end(); ++d ) {
						DT p = T(*k,*d)*E(*d,(int)A[i],(int)B[j])*b(i+1,j+1,*d);
						b(i,j,*k) += p;
						Et(*k,*d) += f(i,j,*k)*p;
					}
					for( vector<int>::const_iterator u = depD(*k,(int)A[i]).begin(); u != depD(*k,(int)A[i]).end(); ++u ) {
						DT p = T(*k,*u)*E(*u,(int)A[i],NOTHING)*b(i+1,j,*u);
						b(i,j,*k) += p;
						Et(*k,*u) += f(i,j,*k)*p;
					}
					for( vector<int>::const_iterator s = depBS(*k).begin(); s != depBS(*k).end(); ++s ) {
						DT p = T(*k,*s)*E(*s,NOTHING,NOTHING)*b(i,j,*s);
						b(i,j,*k) += p;
						Et(*k,*s) += f(i,j,*k)*p;
					}
				}

				// tally emmissions
				for( vector<int>::const_iterator s = depAB((int)A[i-1],NOTHING).begin(); s != depAB((int)A[i-1],NOTHING).end(); ++s )
					Ee(*s,(int)A[i-1],NOTHING) += f(i,j,*s)*b(i,j,*s);
				for( vector<int>::const_iterator s = depAB(NOTHING,(int)B[j-1]).begin(); s != depAB(NOTHING,(int)B[j-1]).end(); ++s )
					Ee(*s,NOTHING,(int)B[j-1]) += f(i,j,*s)*b(i,j,*s);
				for( vector<int>::const_iterator s = depAB((int)A[i-1],(int)B[j-1]).begin(); s != depAB((int)A[i-1],(int)B[j-1]).end(); ++s )
					Ee(*s,(int)A[i-1],(int)B[j-1]) += f(i,j,*s)*b(i,j,*s);
				for( vector<int>::const_iterator s = depAB(NOTHING,NOTHING).begin(); s != depAB(NOTHING,NOTHING).end(); ++s )
					Ee(*s,NOTHING,NOTHING) += f(i,j,*s)*b(i,j,*s);
			}
		}
	}


	void adjust_parameters( blitz::Array<DT,2> &Et, blitz::Array<DT,3> &Ee ) {
		for( int k = 0; k < M; ++k ) {
			// adjust transitions
			DT total = 0.0;
			for( int l = 0; l < M; ++l ) total += Et(k,l);
			for( int l = 0; l < M; ++l ) T(k,l) = Et(k,l)/total;

			// adjust emissions
			total = 0.0;
			for( int x = 0; x <= NOTHING; ++x ) for( int y = 0; y <= NOTHING; ++y ) total += Ee(k,x,y);
			for( int x = 0; x <= NOTHING; ++x ) for( int y = 0; y <= NOTHING; ++y ) {
				//std::cerr << k << '\t' << x << '\t' << y << '\t' << Ee(k,x,y) << '\t' << total << std::endl;
				E(k,x,y) = Ee(k,x,y)/total;
			}
		}
	}

	// Do baum-welch training on D1 and D2 with checkpoints.
	// This will not tally the border back_bd and will tally and compute
	// the left border so long as extend left is true.
	template< typename D1, typename D2 >
	border train_bw_checkpoint( border start_bd, border back_bd, D1 A, D2 B, blitz::Array<DT,2> estimate_T, blitz::Array<DT,3> estimate_E ) const {
		using namespace blitz;
		int n = A.size(), m = B.size();
		Range all = Range::all();

		// store a breakpoint every sqrt(m) positions
		int sqrt_m = (int)floor(sqrt((double)m+1));
		int num_bp = (int)ceil((double)(m+1)/(double)sqrt_m); // due to rounding error we can have slightly more than sqrt_m breakpoints
		data breakpoints(n+1,num_bp,M);

		// working data for computation
		data bpf(n+1,2,M);
		init_forward(start_bd, bpf(all,0,all),A);

		// compute beakpoints and final forward score
		compute_breakpoints( sqrt_m, breakpoints, bpf, A, B );

		// the last region can be odd shaped and needs special initialization
		data f(n+1,sqrt_m+1,M), b(n+1,sqrt_m+1,M); f = 0.0; b = 0.0;

		closed_interval r( sqrt_m*(num_bp-1), m);
		int current_bp = num_bp-1;
		if( r.a == r.b ) { r.a -= sqrt_m; --current_bp; } // if we have a perfect fit then we need to adjust

		init_backward(back_bd,b(all,r.size()-1,all),A);
		while( 1 ) {
			assert( current_bp >= 0 );
			f(all,0,all) = breakpoints(all,current_bp,all);
			D2 adjB( B ); adjB.range.b -= m-r.b; adjB.range.a = adjB.range.b - (r.b - r.a - 1);
			compute_bw( f, b, A, adjB, estimate_T, estimate_E,  r.a != 0 );	// last param check for end of sequence
			if( r.b == m ) tally_backborder( back_bd, f(all,r.size()-1,all), b(all,r.size()-1,all), A, estimate_T, estimate_E );	// tally up the border

			// move region and reset f and b matrices
			r.b = r.a;
			r.a -= sqrt_m;
			if( r.a < 0 ) break;

			b(all,r.size()-1,all) = b(all,0,all);
			b(all,Range(0,r.size()-2),all) = 0.0;
			f = 0.0;

			--current_bp;
		}

		// return border data
		border result_border;
		if( start_bd.extent(firstDim) == 0 ) {
			result_border.resize(1,M);
			result_border(0,all) = b(1,1,all);
		} else {
			result_border.resize( start_bd.shape() );
			//cerr << result_border.shape() << '\t' << b.shape() << endl;
			result_border(all,all) = b(Range(1,result_border.extent(firstDim)),1,all);
		}

		return result_border;
	}

	// Do baum-welch training on D1 and D2 with checkpoints.
	// This will not tally the border back_bd and will tally and compute
	// the left border so long as extend left is true.
	template< typename D1, typename D2 >
	void train_bw( border start_bd, border back_bd, D1 A, D2 B, blitz::Array<DT,2> estimate_T, blitz::Array<DT,3> estimate_E  ) const {
		using namespace blitz;
		int n = A.size(), m = B.size();
		Range all = Range::all();

		// the last region can be odd shaped and needs special initialization
		data f(n+1,m+1,M), b(n+1,m+1,M); f = 0.0; b = 0.0;
		init_forward(start_bd, f(all,0,all),A);
		init_backward(back_bd,b(all,m,all),A);

		compute_bw( f, b, A, B, estimate_T, estimate_E,  false ); // last param check for end of sequence
		//tally_backborder( back_bd, f(all,r.size()-1,all), b(all,r.size()-1,all), A, estimate_T, estimate_E );	// tally up the border
	}



public:
	pair_hmm(int m) : M(m), start_state(0) {
		T.resize(M,M); T = 0.0;
		E.resize(M,ALPHA::SIZE+1,ALPHA::SIZE+1); E = 0.0;
		depS.resize(M);
		depBS.resize(M);
		depAB.resize(ALPHA::SIZE+1,ALPHA::SIZE+1);
		depUL.resize(M,ALPHA::SIZE+1,ALPHA::SIZE+1);
		depL.resize(M,ALPHA::SIZE+1);
		depU.resize(M,ALPHA::SIZE+1);
		depR.resize(M,ALPHA::SIZE+1);
		depD.resize(M,ALPHA::SIZE+1);
		depDR.resize(M,ALPHA::SIZE+1,ALPHA::SIZE+1);
		op_list.resize(NOTHING+1,NOTHING+1);

	}

	inline void set_transition( int i, int j, DT v ) { T(i,j) = v; }
	inline DT get_transition( int i, int j ) const { return T(i,j); }
	inline void set_emmission( int m, int a, int b, DT v ) { E(m,a,b) = v; }
	inline DT get_emmission( int m, int a, int b ) const { return E(m,a,b); }
	void set_start_state( int i ) { start_state = i; }
	void set_end_state( int i ) { end_state = i; }


	struct fa_result {
		int  max_j;
		int  max_i;
		int  best_i;
		int  best_j;
		DT   best_score;
		DT   score;

	};

	// Either left or top should be zero.
	template< class ODDS_MODEL, typename D1, typename D2 >
	fa_result compute_forward_anchor( ODDS_MODEL &odds, D1 A, D2 B ) const {
		using namespace std;
		using namespace blitz;
		int n = (int) A.size(), m = (int) B.size();
		//assert(n == m);

		typename ODDS_MODEL::preprocess_data_ptr odds_data_a = odds.preprocess(A), odds_data_b = odds.preprocess(B);

		Range all = Range::all();

		data f(n+1,2,M);
		border lower(m+1,M), null;
		init_forward(null,f(all,0,all),A);

		fa_result result;
		result.best_j = 0;
		result.best_i = 0;
		result.best_score = 0.0;

		// save first position
		lower(0,all) = f(n,0,all);

		int prev = 0, cur = 1;
		for( int j = 1; j <= m; ++j ) {
			// set everything back to zero
			f(all,cur,all) = 0.0;

			// for i == 0 we can only do L and S transitions
			for( vector<int>::const_iterator k = eval_order.begin(); k != eval_order.end(); ++k ) {
				for( vector<int>::const_iterator l = depL(*k,(int)B[j-1]).begin(); l != depL(*k,(int)B[j-1]).end(); ++l )
					f(0,cur,*k) += f(0,prev,*l)*T(*l,*k)*E(*k,NOTHING,(int)B[j-1]);
				for( vector<int>::const_iterator s = depS(*k).begin(); s != depS(*k).end(); ++s )
					f(0,cur,*k) += f(0,cur,*s)*T(*s,*k)*E(*k,NOTHING,NOTHING);
			}

			// now we cover the remaining blocks
			for( int i = 1; i <= n; ++i ) {
				for( vector<int>::const_iterator k = eval_order.begin(); k != eval_order.end(); ++k ) {
					for( vector<int>::const_iterator l = depL(*k,(int)B[j-1]).begin(); l != depL(*k,(int)B[j-1]).end(); ++l )
						f(i,cur,*k) += f(i,prev,*l)*T(*l,*k)*E(*k,NOTHING,(int)B[j-1]);
					for( vector<int>::const_iterator d = depUL(*k,(int)A[i-1],(int)B[j-1]).begin(); d != depUL(*k,(int)A[i-1],(int)B[j-1]).end(); ++d )
						f(i,cur,*k) += f(i-1,prev,*d)*T(*d,*k)*E(*k,(int)A[i-1],(int)B[j-1]);
					for( vector<int>::const_iterator u = depU(*k,(int)A[i-1]).begin(); u != depU(*k,(int)A[i-1]).end(); ++u )
						f(i,cur,*k) += f(i-1,cur,*u)*T(*u,*k)*E(*k,(int)A[i-1],NOTHING);
					for( vector<int>::const_iterator s = depS(*k).begin(); s != depS(*k).end(); ++s )
						f(i,cur,*k) += f(i,cur,*s)*T(*s,*k)*E(*k,NOTHING,NOTHING);

					ereal score = f(i,cur,*k)/odds.pr(*odds_data_a,0,i-1)/odds.pr(*odds_data_b,0,j-1);
					if( score > result.best_score ) {
						result.best_i = i;
						result.best_j = j;
						result.best_score = score;
					}
				}
			}

			// save bottom position
			lower(j,all) = f(n,cur,all);

			// swap columns
			int tmp = cur; cur = prev; prev = tmp;
		}
		result.max_j = m;
		result.max_i = 0;
		result.score = f(0,prev,0)/odds.pr(*odds_data_b,0,m-1);

		for( int j = 0; j <= m; ++j )
			for( int k = 0; k < M; ++k ) {
				ereal score = lower(j,k)/odds.pr(*odds_data_a,0,n-1)/odds.pr(*odds_data_b,0,j-1);

				if( score > result.score ) {
					result.max_j = j;
					result.max_i = n;
					result.score = score;
				}
			}
		for( int i = 0; i <= n; ++i )
			for( int k = 0; k < M; ++k ) {
				ereal score = f(i,prev,k)/odds.pr(*odds_data_a,0,i-1)/odds.pr(*odds_data_b,0,m-1);

				if( score > result.score ) {
					result.max_j = m;
					result.max_i = i;
					result.score = score;
				}
			}

		return result;
	}

	template<typename SEQA, typename SEQB, typename DP, typename RESULT, typename ODDS, typename ODDSA, typename ODDSB>
	inline void forward_anchor_fill_dp( SEQA &A, SEQB &B, DP &f, RESULT &result, ODDS &odds, ODDSA &odds_data_a, ODDSB &odds_data_b, closed_interval n, closed_interval m ) const {
		using namespace std;
		using namespace blitz;
		using namespace boost;
		Range all = Range::all();

		// set everything in this block to zero
		f(Range(n.a,n.b),Range(m.a,m.b),all) = 0.0;

		ereal score, t;
		for( int j = m.a; j <= m.b; ++j ) {
			for( int i = n.a; i <= n.b; ++i ) {
				for( vector<operation>::const_iterator k = op_list((int)A[i-1],(int)B[j-1]).begin(); k != op_list((int)A[i-1],(int)B[j-1]).end(); ++k ) {
					switch(get<1>(*k)) {
					case MM:
						f(i,j,get<0>(*k)) += f(i-1,j-1,get<2>(*k))*get<3>(*k);
						break;
					case MZ:
						f(i,j,get<0>(*k)) += f(i-1,j,get<2>(*k))*get<3>(*k);
						break;
					case ZM:
						f(i,j,get<0>(*k)) += f(i,j-1,get<2>(*k))*get<3>(*k);
						break;
					case ZZ:
						f(i,j,get<0>(*k)) += f(i,j,get<2>(*k))*get<3>(*k);
					}
				}
				/*
				for( vector<int>::const_iterator k = eval_order.begin(); k != eval_order.end(); ++k ) {
					ereal score = f(i,j,*k)/odds.pr(*odds_data_a,0,i-1)/odds.pr(*odds_data_b,0,j-1);
					if( score > result.best_score ) {
						result.best_i = i;
						result.best_j = j;
						result.best_score = score;
					}
				}*/
				score = f(i,j,end_state)/odds.pr(*odds_data_a,0,i-1)/odds.pr(*odds_data_b,0,j-1);
				if( score > result.best_score ) {
					result.best_i = i;
					result.best_j = j;
					result.best_score = score;
				}
			}
		}
	}

	template<typename SEQA, typename SEQB, typename DP, typename RESULT, typename ODDS, typename ODDSA, typename ODDSB>
	inline void forward_anchor_fill_dpv( SEQA &A, SEQB &B, DP &f, RESULT &result, ODDS &odds, ODDSA &odds_data_a, ODDSB &odds_data_b, closed_interval n, closed_interval m ) const {
		using namespace std;
		using namespace blitz;
		using namespace boost;
		Range all = Range::all();

		// set everything in this block to zero
		//f(Range(n.a,n.b),Range(m.a,m.b),all) = 0.0;

		ereal score = 0.0, t;
		for( int j = m.a; j <= m.b; ++j ) {
			for( int i = n.a; i <= n.b; ++i ) {
				for( vector<operation>::const_iterator k = op_list((int)A[i-1],(int)B[j-1]).begin(); k != op_list((int)A[i-1],(int)B[j-1]).end(); ++k ) {
//					t = f(i+get<1>(*k),j+get<2>(*k),get<3>(*k))*get<4>(*k);
//					if( t > f(i,j,get<0>(*k)) ) f(i,j,get<0>(*k)) = t;
					switch(get<1>(*k)) {
						case MM:
							t = f(i-1,j-1,get<2>(*k))*get<3>(*k);
							break;
						case MZ:
							t = f(i-1,j,get<2>(*k))*get<3>(*k);
							break;
						case ZM:
							t = f(i,j-1,get<2>(*k))*get<3>(*k);
							break;
						case ZZ:
							t = f(i,j,get<2>(*k))*get<3>(*k);
							break;
						case ZERO:
							t = 0.0;
							break;
						default:
							assert(false);
					}
					if( !get<4>(*k) || f(i,j,get<0>(*k)) < t ) f(i,j,get<0>(*k)) = t;

				}
				/*for( vector<int>::const_iterator k = eval_order.begin(); k != eval_order.end(); ++k ) {
					score = f(i,j,*k)/odds.pr(*odds_data_a,0,i-1)/odds.pr(*odds_data_b,0,j-1);
					if( score > result.best_score ) {
						result.best_i = i;
						result.best_j = j;
						result.best_score = score;
					}
				}*/

				score = f(i,j,end_state)/odds.pr(*odds_data_a,0,i-1)/odds.pr(*odds_data_b,0,j-1);
				if( score > result.best_score ) {
					result.best_i = i;
					result.best_j = j;
					result.best_score = score;
				}

			}
		}
	}


	template<typename DP, typename RESULT, typename ODDS, typename ODDSA, typename ODDSB>
	inline void forward_anchor_eval_max( DP &f, RESULT &result, ODDS &odds, ODDSA &odds_data_a, ODDSB &odds_data_b, int n, int m ) const {
		using namespace std;
		using namespace blitz;
		Range all = Range::all();

		result.max_j = m;
		result.max_i = 0;
		result.score = f(0,m,0)/odds.pr(*odds_data_b,0,m-1);

		for( int j = 0; j <= m; ++j ) {
			//for( int k = 0; k < M; ++k ) {
				ereal score = f(n,j,end_state)/odds.pr(*odds_data_a,0,n-1)/odds.pr(*odds_data_b,0,j-1);
				if( score > result.score ) {
					result.max_j = j;
					result.max_i = n;
					result.score = score;
				}
			//}
		}
		for( int i = 0; i <= n; ++i ) {
			//for( int k = 0; k < M; ++k ) {
				ereal score = f(i,m,end_state)/odds.pr(*odds_data_a,0,i-1)/odds.pr(*odds_data_b,0,m-1);
					if( score > result.score ) {
					result.max_j = m;
					result.max_i = i;
					result.score = score;
				}
			//}
		}
	}

	// Either left or top should be zero.
	template< typename DP, class ODDS_MODEL, typename D1, typename D2 >
	fa_result compute_forward_anchor_b( DP &f, ODDS_MODEL &odds, D1 A, D2 B, int step = 50, double drop_threshold = -50, double ext_threshold = 10.0, bool fast = true ) const {
		using namespace std;
		using namespace blitz;
		int n = (int) A.size(), m = (int) B.size();
		//assert(n == m);

		typename ODDS_MODEL::preprocess_data_ptr odds_data_a = odds.preprocess(A), odds_data_b = odds.preprocess(B);

		Range all = Range::all();

//		data f(n+1,m+1,M);
		border null;
		init_forward(null,f(all,0,all),A);
		init_forward(null,f(0,all,all),B);

		fa_result result;
		result.best_j = 0;
		result.best_i = 0;
		result.best_score = 0.0;

		int cn = min(n,step), cm = min(m,step);

		// fill first block
		closed_interval inter_n(1,cn), inter_m(1,cm);


		if(fast) forward_anchor_fill_dpv( A, B, f, result, odds, odds_data_a, odds_data_b, inter_n, inter_m ) ;
		else     forward_anchor_fill_dp( A, B, f, result, odds, odds_data_a, odds_data_b, inter_n, inter_m );


		while(1) {
			forward_anchor_eval_max( f, result, odds, odds_data_a, odds_data_b, cn, cm );

			if( result.score.as_base() >= ext_threshold || result.score.as_base() <= drop_threshold ) return result;

			// can't expand horizon
			if( cn == n && cm == m ) break;

			// check and expand in n dimension
			if( cn < n ) {
				closed_interval new_n( cn+1, min(n,cn+step) );
				closed_interval new_m( 1, cm );
				if(fast) forward_anchor_fill_dpv( A, B, f, result, odds, odds_data_a, odds_data_b, new_n, new_m );
				else     forward_anchor_fill_dp( A, B, f, result, odds, odds_data_a, odds_data_b, new_n, new_m );
			}

			// check expand in dimension m
			if( cm < m ) {
				closed_interval new_m( cm+1, min(m,cm+step) );
				closed_interval new_n( 1, min(n,cn+step) );
				if(fast) forward_anchor_fill_dpv( A, B, f, result, odds, odds_data_a, odds_data_b, new_n, new_m );
				else     forward_anchor_fill_dp( A, B, f, result, odds, odds_data_a, odds_data_b, new_n, new_m );
			}

			cn = min(n,cn+step);
			cm = min(m,cm+step);
		}


		return result;
	}



	typedef blitz::Array< DT, 2 > dprow;
	dprow rowA, rowB;
	struct row_metrics {
		closed_interval valid_region;
		DT best_score;
		int score_col;
		row_metrics() : valid_region(-1,-1), best_score(0.0), score_col(-1) {}
	};


	template< class DPROW >
	inline void dp_forward_init_cell( int i, DPROW &crow ) {
		using namespace boost;
		using namespace std;
		using namespace blitz;
		crow(i,Range::all()) = 0.0;
		crow(i,start_state) = 1.0;

		for( vector<int>::const_iterator k = eval_order.begin(); k != eval_order.end(); ++k )
			for( vector<int>::const_iterator s = depS(*k).begin(); s != depS(*k).end(); ++s )
				crow(i,*k) += crow(i,*s)*T(*s,*k)*E(*k,NOTHING,NOTHING);
	}

	template< class DPROW, typename D2, typename O2 >
	inline row_metrics dp_forward_row(
			closed_interval fill_region,
			DPROW &crow,
			DPROW &prow,
			std::vector<op_vector> &op_row,
			D2 &A,
			O2 &oddsA,
			DT oddsB,
			double drop_threshold,
			bool first_row = false ) {

		using namespace blitz;
		using namespace btl;
		using namespace boost;

		Range all = Range::all();

		row_metrics rmetrics;
		rmetrics.valid_region.a = fill_region.a;

		DT score, t;

		for( int i = fill_region.a; i <= fill_region.b; ++i ) {
			assert( (int)A[i-1] <= 5 );
			for( op_vector::const_iterator k = op_row[(int)A[i-1]].begin(); k != op_row[(int)A[i-1]].end(); ++k ) {
				switch(get<1>(*k)) {
					case MM:
						//cerr << "MM" << "\t" << get<2>(*k) << "\t" << prow(i-1,get<2>(*k)) << "\t" << get<3>(*k) << endl;
						t = prow(i-1,get<2>(*k))*get<3>(*k);
						break;
					case MZ:
						//cerr << "MZ" << "\t" << get<2>(*k) << "\t" << crow(i-1,get<2>(*k)) << "\t" << get<3>(*k) << endl;
						t = crow(i-1,get<2>(*k))*get<3>(*k);
						break;
					case ZM:
						//cerr << "ZM" << "\t" << get<2>(*k) << "\t" << prow(i,get<2>(*k)) << "\t" << get<3>(*k) << "\t" << &prow << endl;
						t = prow(i,get<2>(*k))*get<3>(*k);
						break;
					case ZZ:
						//cerr << "ZZ" << "\t"<< get<2>(*k) << "\t" << crow(i,get<2>(*k)) << "\t" << get<3>(*k) << endl;
						t = crow(i,get<2>(*k))*get<3>(*k);
						break;
					case ZERO:
						t = 0.0;
						break;
					default:
						assert(false);
				}
#ifdef VITERBI_EXTENSIONS
				if( !get<4>(*k)  ) crow(i,get<0>(*k)) = t; else crow(i,get<0>(*k)) = max(crow(i,get<0>(*k)),t);
#else
				if( !get<4>(*k)  ) crow(i,get<0>(*k)) = t; else crow(i,get<0>(*k)) += t;
#endif

				//cerr << "ASN " << i << "\t" << get<0>(*k) << "\t" << crow(i,get<0>(*k)).as_base() << "\t" << (int)&crow << endl;
			}

			score = crow(i,end_state)/oddsA.pr( A.global_coord( closed_interval(0,i-1) ) )/oddsB;
			//cerr << i << "\t" << score.as_base() << "\t" << crow(i,end_state)
			//<< "\t" <<  oddsA.pr(A.global_coord( closed_interval(0,i-1) )) << endl;

			if( score > rmetrics.best_score ) {
				rmetrics.best_score = score;
				rmetrics.score_col = i;
			}
			if( score.as_base() <= drop_threshold && i == rmetrics.valid_region.a ) ++rmetrics.valid_region.a;
			else if( score.as_base() > drop_threshold ) rmetrics.valid_region.b = i;
			if( first_row  && score.as_base() <= drop_threshold ) break;
		}
		return rmetrics;
	}

	/** New forward extension algorithm inspired by blast's x-drop.
	 *
	 *	The two input sequences are assumed to be appropriately left bounded
	 *	sequence regions. Reverse extensions can be done using the reverse
	 *	sequence region data type.
	 */
	template< typename O1, typename O2, typename D1, typename D2, typename OutputIterator, typename Bounds >
	std::pair<int,DT> forward_extension( D1 &A, O1 &oddsA, D2 &B, O2 &oddsB, double drop_threshold, OutputIterator out, Bounds &bounds ) {
		using namespace std;
		using namespace boost;
		using namespace blitz;
		Range all = Range::all();

		// We keep the row dp data arrays allocated for performance. We allocate the data
		// lazily here. We could allocate only enough for the range, but then we'd be
		// reallocating often for no reason.
		if( rowA.extent(firstDim) < (int)A.data.size()+1 ) {
			rowA.resize(A.data.size()+1,M);
			rowB.resize(A.data.size()+1,M);
		}

#ifdef DEBUG_OUTPUT
		unsigned long int fill_area = 0, tail_area = 0;
#endif
		bool current_row = true;
		closed_interval fill_region(1,A.data.size()), tail_region(1,A.data.size()), clipped_fill_region(-1,-1), explored_region(-1,-1);

		// Initialize the first row.
		dp_forward_init_cell( 0, rowA );
		row_metrics rm = dp_forward_row( fill_region, rowA, rowB, operations[NOTHING], A.data, oddsA, 1.0,  drop_threshold, true );
		fill_region = rm.valid_region;
		DT best_score = rm.best_score, coddsB = 1.0;
		int j = 1, max_row = 0;

		dprow *prev, *cur;
		for( ; j <= (int)B.data.size() && fill_region.a <= fill_region.b; ++j ) {
			// Clip fill region.
			clipped_fill_region = bounds.clip( A.data.global_coord(closed_interval(fill_region.a-1,fill_region.b-1))  );
			clipped_fill_region = A.data.local_coord( clipped_fill_region );
			clipped_fill_region.a++; clipped_fill_region.b++;
			explored_region = clipped_fill_region;

			// Break if we are clipped out of a valid fill region.
			if( clipped_fill_region.a > clipped_fill_region.b ) break;

			// Get odds for row.
			coddsB = oddsB.pr( B.data.global_coord( closed_interval(0,j-1) ) );

			// Assign correct row pointers.
			if( current_row ) {
				cur = &rowB;
				prev = &rowA;
			} else {
				cur = &rowA;
				prev = &rowB;
			}

			// Set left-most cell to zero.
			(*cur)(clipped_fill_region.a-1,all) = 0.0;

			double drop = max( 0.0 + drop_threshold, best_score.as_base() + drop_threshold );
			rm = dp_forward_row( clipped_fill_region, *cur, *prev, operations[(int)B.data[j-1]], A.data, oddsA, coddsB,  drop, false );
			//cerr << "pre " << rm.valid_region << "\t" << rm.best_score.as_base() << "\t" << rm.score_col << endl;

			// We may need to extend past the previous row boundary.
			if( rm.valid_region.b == fill_region.b ) {
				tail_region.a = fill_region.b+1;
				clipped_fill_region = bounds.clip( A.data.global_coord(closed_interval(tail_region.a-1,tail_region.b-1))  );
				clipped_fill_region = A.data.local_coord( clipped_fill_region );
				clipped_fill_region.a++; clipped_fill_region.b++;

				// Make sure tail is not clipped.
				if( clipped_fill_region.a <= clipped_fill_region.b  && clipped_fill_region.a == tail_region.a ) {
					row_metrics rme = dp_forward_row( clipped_fill_region, *cur, *prev, operations[NOTHING], A.data, oddsA, coddsB, drop, true );

					// If end of valid region was right on the border we can get strange
					// results here. So we check that the results are consistent with a
					// valid tail.
					if( rme.valid_region.a == tail_region.a && rme.valid_region.b > tail_region.a ) {
						explored_region.b = rme.valid_region.b;
						rm.valid_region.b = rme.valid_region.b;
						if( rm.best_score < rme.best_score ) {
							rm.best_score = rme.best_score;
							rm.score_col = rme.score_col;
						}
					}
				}
			}
			//cerr << j << "\t" << rm.valid_region << "\t" << rm.best_score.as_base() << "\t" << rm.score_col << endl;

			// Save extension data.
			--explored_region.a; --explored_region.b;
			out = extension_store::extension::row( A.data.global_coord(explored_region), A.data.global_coord(rm.score_col - 1) );
			++out;

			// Set new boundaries and adjust max score.
			fill_region = rm.valid_region;

#ifdef DEBUG_OUTPUT
			tail_area += rm.valid_region.size();
#endif

			if( rm.best_score > best_score ) {
				best_score = rm.best_score;
				max_row = j;
#ifdef DEBUG_OUTPUT
				fill_area += tail_area;
				tail_area = 0;
#endif
			}

			// Step bounds.
			bounds.step( A.data.global_coord(closed_interval(fill_region.a-1,fill_region.b-1)) );

			// Flip current_row.
			current_row = current_row ^ true;
		}
#ifdef DEBUG_OUTPUT
		cerr << "fill_area " << fill_area << endl;
		cerr << "tail_area " << tail_area << endl;
#endif
		return make_pair(max_row - 1, best_score);
	}

	template< class DPROW, typename D1, typename O1, typename D2, typename O2 >
	inline row_metrics dp_ungapped_forward_row( DPROW &row, D1 &A, O1 &oddsA, D2 &B, O2 &oddsB, double drop_threshold )	{
		using namespace blitz;
		using namespace btl;
		using namespace boost;

		Range all = Range::all();

		row_metrics rmetrics; rmetrics.valid_region.a = 1;

		DT score, t;

		for( int i = 1, j = 1; i <= (int)A.size() && j <= (int)B.size(); ++i, ++j ) {
			for( op_vector::const_iterator k = operations[(int)B[j-1]][(int)A[i-1]].begin(); k != operations[(int)B[j-1]][(int)A[i-1]].end(); ++k ) {
				switch(get<1>(*k)) {
					case MM:
						//cerr << "MM" << "\t" << get<2>(*k) << "\t" << row(i-1,get<2>(*k)).as_base() << "\t" << get<3>(*k) << endl;
						t = row(i-1,get<2>(*k))*get<3>(*k);
						break;
					case MZ:
						t = 0.0;
						break;
					case ZM:
						t = 0.0;
						break;
					case ZZ:
						//cerr << "ZZ" << "\t"<< get<2>(*k) << "\t" << row(i,get<2>(*k)) << "\t" << get<3>(*k) << endl;
						t = row(i,get<2>(*k))*get<3>(*k);
						break;
					case ZERO:
						t = 0.0;
						break;
					default:
						assert(false);
				}
				//cerr << "ASN " << i << "\t" << get<0>(*k) << "\t" << row(i,get<0>(*k)).as_base() << "\t" << t.as_base() << "\t" << get<4>(*k)<< endl;
				if( !get<4>(*k)  ) row(i,get<0>(*k)) = t; else row(i,get<0>(*k)) += t;
			}
			score = row(i,end_state)/oddsA.pr( A.global_coord( closed_interval(0,i-1) ) )/oddsB.pr( B.global_coord( closed_interval(0,j-1)) );
			//cerr << i << "\t" << score.as_base() << "\t" << row(i,end_state)
			//j<< "\t" <<  oddsA.pr(A.global_coord( closed_interval(0,i-1) )).as_base()
			//<< '\t' << oddsB.pr( B.global_coord( closed_interval(0,j-1)) ).as_base() << endl;
			//cerr << rmetrics.best_score << "\t" << score << endl;
			if( score > rmetrics.best_score ) {
				rmetrics.best_score = score;
				rmetrics.score_col = i;
			}
			if( score.as_base() > drop_threshold ) rmetrics.valid_region.b = i;
			else break;

		}
		return rmetrics;
	}

	/** Ungapped forward extension.
	 */
	template< typename O1, typename O2, typename D1, typename D2 >
	std::pair<DT, int> ungapped_forward_extension( D1 &A, O1 &oddsA, D2 &B, O2 &oddsB, double drop_threshold ) {
		using namespace std;
		using namespace boost;
		using namespace blitz;
		Range all = Range::all();

		// We keep the row dp data arrays allocated for performance. We allocate the data
		// lazily here. We could allocate only enough for the range, but then we'd be
		// reallocating often for no reason.
		unsigned int buf_size = max(A.data.size(),B.data.size()) + 1;	
		if( rowA.extent(firstDim) < buf_size ) {
			rowA.resize(buf_size,M);
			rowB.resize(buf_size,M);
		}

		// Initialize the first row.
		dp_forward_init_cell( 0, rowA );
		row_metrics rm = dp_ungapped_forward_row( rowA, A.data, oddsA, B.data, oddsB, drop_threshold );
		return make_pair( rm.best_score, rm.score_col );
	}

	// Do baum-welch training on regions of D1 and D2.
	// Expects two iterators, start and end, pointing to region pairs.
	template< typename D1, typename D2, typename ITERATOR >
	DT constrained_train_bw( ITERATOR start, ITERATOR end, int interface_size = 4, bool status = false ) {
		using namespace std;
		using namespace blitz;
		Range all = Range::all();

		// first compute all borders
		std::vector<border> forward_borders;

		ITERATOR next = start, it = start; ++next;
		D1 firstRa = *(it->first);
		firstRa.data.range.b = firstRa.data.range.a + interface_size*2 - 2;
		border first_border(interface_size*2,M), nullbd;
		init_forward( border(), first_border, firstRa.data);
		forward_borders.push_back( first_border );

		data f;

		while(1) {
			// first adjust region to align based on interface
			D1 rega = *(it->first);	D2 regb = *(it->second);

			// if not at end, add lower interface
			if( next != end ) rega.data.range.b += interface_size - 1;

			// if not at start, add upper interface. Previous region included "interface" bases from a.
			if( it != start ) rega.data.range.a -= interface_size;

			// allocate dp matricies and fill
			int n = (int) rega.data.size();
			f.resize(n+1,2,M); f = 0.0;

			init_forward(forward_borders.back(),f(all,0,all),rega.data);
			//if( status ) std::cerr << "Doing forward computation " << rega.data.range << " vs " << regb.data.range << std::endl;

			// compute breakpoints and final forward score
			int col = compute_breakpoints( 0, data(), f, rega.data, regb.data );

			// copy border data back to idata
			forward_borders.push_back( f(Range(n-interface_size*2+1,n),col,all).copy() );
			++next;
			++it;
			if( it == end ) break;
		}

		// save final score
		DT final_sc = forward_borders.back()(interface_size*2-1,end_state);
		forward_borders.pop_back();

		// now go backwards doing bw
		blitz::Array<DT,2> estimate_T(M,M);	estimate_T = 0.0;
		blitz::Array<DT,3> estimate_E(M,NOTHING+1,NOTHING+1);	estimate_E = 0.0;

		it = end; --it; next = end;
		D1 lastRa = *(it->first);
		lastRa.data.range.a = lastRa.data.range.b - interface_size*2 + 1;
		border backward_border_init(interface_size*2+1,M);
		init_backward( border(), backward_border_init, lastRa.data );
		border backward_border = backward_border_init(Range(1,interface_size*2),all);

		while(1) {
			// setup region
			D1 rega = *(it->first);	D2 regb = *(it->second);

			// if not at front, add interface at front
			if( it != start ) rega.data.range.a -= interface_size;

			// if not at end, add interface at end
			if( next != end ) {
				rega.data.range.b += interface_size;
				regb.data.range.b += 1;
			}

			//if( status ) std::cerr << "Training " << rega.data.range << " vs " << regb.data.range << std::endl;
			//cerr << forward_borders.back() << endl;
			//cerr << backward_border << endl;
			backward_border = train_bw_checkpoint( forward_borders.back(), backward_border, rega.data, regb.data, estimate_T, estimate_E );

			// move to next forward border
			forward_borders.pop_back();

			if( it == start ) break;

			--it;
			--next;
		}

		// adjust parameters with learned values
		adjust_parameters( estimate_T, estimate_E );
		return final_sc;
	}


	template< class DPROW, class PATHROW, typename D2>
		inline void vt_row(
				closed_interval fill_region,
				DPROW &crow,
				DPROW &prow,
				PATHROW &cpath,
				PATHROW &ppath,
				std::vector<op_vector> &op_row,
				D2 &A,
				int B
			) {

			using namespace blitz;
			using namespace btl;
			using namespace boost;

			Range all = Range::all();

			DT score, t;
			for( int i = fill_region.a; i <= fill_region.b; ++i ) {
				for( op_vector::const_iterator k = op_row[(int)A[i-1]].begin(); k != op_row[(int)A[i-1]].end(); ++k ) {
					switch(get<1>(*k)) {
						case MM:
							//cerr << "MM" << "\t" << get<2>(*k) << "\t" << prow(i-1,get<2>(*k)) << "\t" << get<3>(*k) << endl;
							t = prow(i-1,get<2>(*k))*get<3>(*k);
							if( !get<4>(*k) || t > crow(i,get<0>(*k)) ) {
								crow(i,get<0>(*k)) = t;
								cpath(i,get<0>(*k)).set( ppath(i-1,get<2>(*k)), path_column((int)A[i-1],B,get<0>(*k)) );
							}
							break;
						case MZ:
							//cerr << "MZ" << "\t" << get<2>(*k) << "\t" << crow(i-1,get<2>(*k)) << "\t" << get<3>(*k) << endl;
							t = crow(i-1,get<2>(*k))*get<3>(*k);
							if( !get<4>(*k) || t > crow(i,get<0>(*k)) ) {
								crow(i,get<0>(*k)) = t;
								cpath(i,get<0>(*k)).set( cpath(i-1,get<2>(*k)), path_column((int)A[i-1],NOTHING,get<0>(*k)) );
							}

							break;
						case ZM:
							//cerr << "ZM" << "\t" << get<2>(*k) << "\t" << prow(i,get<2>(*k)) << "\t" << get<3>(*k) << "\t" << &prow << endl;
							t = prow(i,get<2>(*k))*get<3>(*k);
							if( !get<4>(*k) || t > crow(i,get<0>(*k)) ) {
								crow(i,get<0>(*k)) = t;
								cpath(i,get<0>(*k)).set( ppath(i,get<2>(*k)), path_column(NOTHING,B,get<0>(*k)) );
							}
							break;
						case ZZ:
							//cerr << "ZZ" << "\t"<< get<2>(*k) << "\t" << crow(i,get<2>(*k)) << "\t" << get<3>(*k) << endl;
							t = crow(i,get<2>(*k))*get<3>(*k);
							if( !get<4>(*k) || t > crow(i,get<0>(*k)) ) {
								crow(i,get<0>(*k)) = t;
								cpath(i,get<0>(*k)).set( cpath(i,get<2>(*k)), path_column(NOTHING,NOTHING,get<0>(*k)) );
							}
							break;
						case ZERO:
							t = 0.0;
							break;
						default:
							assert(false);
					}
				//cerr << "ASN " << i << "\t" << get<0>(*k) << "\t" << crow(i,get<0>(*k)).as_base() << "\t" << (int)&crow << endl;
				}
			}
		}

	// We store our global training information here.
	DT final_pr;
	blitz::Array<DT,2> eT;
	blitz::Array<DT,3> eE;

	void start_training() {
		final_pr = 1.0;
		eT.resize(M,M); eT = 0.0;
		eE.resize(M,NOTHING+1,NOTHING+1); eE = 0.0;
	}

	/** Do Baum-Welch training on over a set of regions between D1 and D2.
	 *
	 *  The iterators start and end define a set of segment pairs sequences. Each
	 *  segment pair sequence is defined by
	 *  We collect data over each sequence of segment pairs and adjust parameters
	 *  for the whole set of sequences.
	 */
	// Expects two iterators, start and end, pointing to region pairs.
	template< typename D1, typename D2, typename ITERATOR>
	DT constrained_train_bw_mt( int nthreads, ITERATOR start, ITERATOR end, int interface_size = 4, bool status = false ) {
		using namespace std;
		using namespace blitz;
		Range all = Range::all();
//		if( nthreads < 2 ) return constrained_train_bw<D1,D2>( start, end, interface_size, status );

		// package data for workers
		typedef typename pair_hmm_bw_worker<ALPHA,DT,D1,D2,ITERATOR>::worker_data worker_data;
		typename pair_hmm_bw_worker<ALPHA,DT,D1,D2,ITERATOR>::work_queue wqueue;

		DT pr = 1.0;

		// first compute all borders
		std::vector<border> forward_borders;
		std::queue<border> backward_borders;
		pair_hmm_forward_worker<ALPHA,DT,D1,D2,ITERATOR> fw_worker( *this, start, end, interface_size, forward_borders );
		pair_hmm_backward_worker<ALPHA,DT,D1,D2,ITERATOR> bw_worker( *this, start, end, interface_size, backward_borders );
		{
			boost::thread_group workers;
			workers.create_thread(fw_worker);
			workers.create_thread(bw_worker);
			workers.join_all();
		}

		// save final score
		pr = forward_borders.back()(1,end_state);
		final_pr *= pr;
		forward_borders.pop_back();

		//ITERATOR lastit = i->second;
		//--lastit;
		//cerr << "Pr-sub: " << pr << '\t' << i->first->first->data.range.a << '\t';
		//cerr << lastit->first->data.range.b << '\t';
		//cerr << i->first->second->data.range.a << '\t' << lastit->second->data.range.b << endl;


		ITERATOR it = end, next = end; --it;
		while(1) {
			// setup region
			D1 rega = *(it->first);	D2 regb = *(it->second);
			//cerr << it->first->data.range << '\t';
			//cerr << it->second->data.range << '\n';


			// if not at front, add interface at front
			if( it != start ) rega.data.range.a -= interface_size;

			// if not at end, add interface at end
			if( next != end ) rega.data.range.b += interface_size;

			//cerr << rega.data.range << '\t' << regb.data.range << endl;
			//cerr << forward_borders.size() << '\t' <<backward_borders.size() << endl;

			assert( forward_borders.size() > 0 );
			assert( backward_borders.size() > 0 );

			wqueue.push_back( worker_data( forward_borders.back(), backward_borders.front(), rega, regb ) );

			// move to next forward border
			forward_borders.pop_back();
			backward_borders.pop();

			if( it == start ) break;

			--it;
			--next;
		}
		wqueue.shutdown_when_empty();

		// Create workers.
		blitz::Array<DT,3> estimate_T(nthreads,M,M); estimate_T = 0.0;
		blitz::Array<DT,4> estimate_E(nthreads,M,NOTHING+1,NOTHING+1);	estimate_E = 0.0;
		{
			boost::thread_group workers;
			for( int i = 0; i < nthreads; ++i ) {
				pair_hmm_bw_worker<ALPHA,DT,D1,D2,ITERATOR> worker( *this, wqueue, estimate_T(i,all,all), estimate_E(i,all,all,all) );
				workers.create_thread( worker );
			}
			workers.join_all();
		}

		for( int i = 0; i < nthreads; ++i ) {
			for( int j = 0; j < M; ++j ) {
				for( int k = 0; k < M; ++k ) {
					eT(j,k) += estimate_T(i,j,k) / pr;
				}
				for( int x = 0; x <= NOTHING; ++x ) {
					for( int y = 0; y <= NOTHING; ++y ) {
						eE(j,x,y) += estimate_E(i,j,x,y) / pr;
					}
				}
			}
		}
		return pr;
	}

	DT end_training() {
		// adjust parameters with learned values
		adjust_parameters( eT, eE );
		return final_pr;
	}


	// do baum-welch training on D1 and D2
	template< typename D1, typename D2 >
	DT train_bw( D1 A, D2 B ) {
		using namespace std;
		using namespace blitz;

		Array<DT,2> estimate_T(M,M);	estimate_T = 0.0;
		Array<DT,3> estimate_E(M,NOTHING+1,NOTHING+1);	estimate_E = 0.0;

		border result = train_bw_checkpoint( border(), border(), A, B, estimate_T, estimate_E );

		// adjust parameters with learned values
		adjust_parameters( estimate_T, estimate_E );

		assert( result.extent(firstDim) == 1 );
		return result(0,end_state);
	}

	// Output all probabilities.
	std::ostream &write_parameters( std::ostream &out ) {
		// output transition probabilities
		for( int i = 0; i < M; ++i ) {
			for( int j = 0; j < M; ++j ) {
				out << (double) T(i,j) << " ";
			}
			out << std::endl;
		}

		// output emmision probabilities
		for( int i = 0; i < M; ++i ) {
			for( int j = 0; j <= NOTHING; ++j ) {
				for( int k = 0; k <= NOTHING; ++k ) {
					out << (double) E(i,j,k) << " ";
				}
				out << std::endl;
			}
			out << std::endl;
		}
		return out;
	}

	std::istream &read_parameters( std::istream &in ) {
		return in;
	}

	// determine state evaluation order and dependancies to ignore zero'd entries
	void optimize() {
		using namespace std;
		using namespace boost;

		// clear and resize existing lists
		eval_order.clear();
		eval_orderB.clear();
		depS = std::vector<int>();
		depBS = std::vector<int>();
		depL = std::vector<int>();
		depR = std::vector<int>();
		depU = std::vector<int>();
		depD = std::vector<int>();
		depUL = std::vector<int>();
		depDR = std::vector<int>();
		depAB = std::vector<int>();

		// build dependancy lists
		for( int j = 0; j < M; ++j ) {
			for( int i = 0; i < M; ++i ) {
				if( T(i,j) > 0.0 && E(j,NOTHING,NOTHING) > 0.0 ) {
					depS(j).push_back(i);
					depBS(i).push_back(j);
				}

				for( int x = 0; x < NOTHING; ++x ) {
					if( T(i,j) > 0.0 && E(j,NOTHING,x) > 0.0 ) {
						depL(j,x).push_back(i);
						depR(i,x).push_back(j);
					}
					if( T(i,j) > 0.0 && E(j,x,NOTHING) > 0.0 ) {
						depU(j,x).push_back(i);
						depD(i,x).push_back(j);
					}
					for( int y = 0; y < NOTHING; ++y ) {
						if( T(i,j) > 0.0 && E(j,x,y) > 0.0 ) {
							depUL(j,x,y).push_back(i);
							depDR(i,x,y).push_back(j);
						}
					}
				}
			}
		}

		for( int j = 0; j < M; ++j )
			for( int x = 0; x <= NOTHING; ++x )
				for( int y = 0; y <= NOTHING; ++y )
					if( E(j,x,y) > 0.0 ) depAB(x,y).push_back(j);


		// catagorize states in to silent and non-silent states
		set<int> silent;
		for( int i = 0; i < M; ++i ) if( depS(i).size() > 0.0 ) silent.insert(i); else eval_order.push_back(i);

		// work out proper eval order
		while( silent.size() > 0 ) {
			unsigned int old_size = silent.size();

			for( set<int>::iterator i = silent.begin(); i != silent.end(); ++i ) {
				bool valid = true;
				// if this state depends on a silent state not yet processed, we can't add it yet
				for( vector<int>::const_iterator j = depS(*i).begin(); j != depS(*i).end(); ++j )
					if( silent.find(*j) != silent.end() ) {
						valid = false;
						break;
					}

				// skip to the next element if this one is not valid
				if( !valid ) continue;

				// otherwise add it to the order, remove it from the set, and search again
				eval_order.push_back(*i);
				silent.erase(i);
				break;
			}

			// no progress was made, silent cycle in the model
			if( old_size == silent.size() ) {
				std::cerr << "pair_HMM: silent cycle in model" << std::endl;
				for( set<int>::iterator i = silent.begin(); i != silent.end(); ++i ) {
					std::cerr << "silent state " << *i << " depends on" << std::endl;
					// if this state depends on a silent state not yet processed, we can't add it yet
					for( vector<int>::const_iterator j = depS(*i).begin(); j != depS(*i).end(); ++j )
						std::cerr << *j << std::endl;
				}
				exit(1);
			}
		}

		// catagorize states in to silent and non-silent states
		silent.clear();
		for( int i = 0; i < M; ++i ) if( depBS(i).size() > 0 ) silent.insert(i); else eval_orderB.push_back(i);

		// work out proper eval order
		while( silent.size() > 0 ) {
			unsigned int old_size = silent.size();

			for( set<int>::iterator i = silent.begin(); i != silent.end(); ++i ) {
				bool valid = true;
				// if this state depends on a silent state not yet processed, we can't add it yet
				for( vector<int>::const_iterator j = depBS(*i).begin(); j != depBS(*i).end(); ++j )
					if( silent.find(*j) != silent.end() ) {
						valid = false;
						break;
					}

				// skip to the next element if this one is not valid
				if( !valid ) continue;

				// otherwise add it to the order, remove it from the set, and search again
				eval_orderB.push_back(*i);
				silent.erase(i);
				break;
			}

			// no progress was made, silent cycle in the model
			if( old_size == silent.size() ) {
				std::cerr << "pair_HMM: silent cycle in model" << std::endl;
				for( set<int>::iterator i = silent.begin(); i != silent.end(); ++i ) {
					std::cerr << "-- " << *i << std::endl;
					// if this state depends on a silent state not yet processed, we can't add it yet
					for( vector<int>::const_iterator j = depBS(*i).begin(); j != depBS(*i).end(); ++j )
						std::cerr << *j << std::endl;
				}
				exit(1);
			}
		}

		// This is in the format j, i, assuming we iterate along i.
		operations.clear();
		operations.resize( NOTHING+1, vector<op_vector>( NOTHING+1 ) );

		blitz::Array<bool,1> seen(M), seen2(M);
		op_list = std::vector<operation>();
		for( int x = 0; x <= NOTHING; ++x ) {
			for( int y = 0; y <= NOTHING; ++y ) {
				seen = false;
				std::vector<operation> &opvec = op_list(x,y);
				op_vector &opvec2 = operations[y][x];
				for( vector<int>::const_iterator k = eval_order.begin(); k != eval_order.end(); ++k ) {
					for( vector<int>::const_iterator l = depL(*k,y).begin(); l != depL(*k,y).end(); ++l ) {
						opvec.push_back( operation(*k,ZM,*l, T(*l,*k)*E(*k,NOTHING,y), seen(*k) ) );
						opvec2.push_back( operation(*k,ZM,*l, T(*l,*k)*E(*k,NOTHING,y), seen(*k) ) );
						seen(*k) = true;
					}
					for( vector<int>::const_iterator d = depUL(*k,x,y).begin(); d != depUL(*k,x,y).end(); ++d ) {
						opvec.push_back( operation(*k,MM,*d, T(*d,*k)*E(*k,x,y), seen(*k) ) );
						opvec2.push_back( operation(*k,MM,*d, T(*d,*k)*E(*k,x,y), seen(*k) ) );
						seen(*k) = true;
					}
					for( vector<int>::const_iterator u = depU(*k,x).begin(); u != depU(*k,x).end(); ++u ) {
						opvec.push_back( operation(*k,MZ,*u, T(*u,*k)*E(*k,x,NOTHING), seen(*k) ) );
						opvec2.push_back( operation(*k,MZ,*u, T(*u,*k)*E(*k,x,NOTHING), seen(*k) ) );
						seen(*k) = true;
					}
					for( vector<int>::const_iterator s = depS(*k).begin(); s != depS(*k).end(); ++s ) {
						opvec.push_back( operation(*k,ZZ,*s,T(*s,*k)*E(*k,NOTHING,NOTHING), seen(*k) ) );
						opvec2.push_back( operation(*k,ZZ,*s,T(*s,*k)*E(*k,NOTHING,NOTHING), seen(*k) ) );
						seen(*k) = true;
					}
					// If *k is never assigned, just set it to zero.
					if( !seen( *k ) ) {
						opvec.push_back( operation(*k,ZERO,*k,0.0, seen(*k) ) );
						opvec2.push_back( operation(*k,ZERO,*k,0.0, seen(*k) ) );
					}
				}
			}
		}

		/*
		cerr << "optimize" << endl;
		int x = 0, y = NOTHING;
		for( unsigned int i = 0;  i < op_list(x,y).size();  ++i ) {
			cerr << get<0>(op_list(x,y)[i]) << "\t";
			cerr << (int) get<1>(op_list(x,y)[i]) << "\t";
			cerr << get<2>(op_list(x,y)[i]) << "\t";
			cerr << get<3>(op_list(x,y)[i]) << "\t";
			cerr << get<4>(op_list(x,y)[i]) << endl;
		}
		*/


	}

};

#endif
