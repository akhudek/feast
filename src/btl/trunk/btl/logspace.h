/*
 *  logspace.h
 *
 *  Created by Alexander K. Hudek on 2008-08-05.
 *
 *	 Stores floating point values in log space. Does not support subtraction or negative numbers.
 */
#ifndef BTL_LOGSPACE_H__
#define BTL_LOGSPACE_H__

#include <cmath>
#include <limits>
#include <iostream>

namespace btl {

	// static lookup tables
	template< typename FPT >
	class logspace {
	public:
		static FPT logadd[16777216];
		static unsigned int const bits;
		static unsigned int const size;
		static FPT const logadd_cutoff;
		static FPT const scale;

		// compute tables
		static void init();

	private:
		double value;

	public:

		// Constructors for various types.
		logspace() {}
		logspace( float v ) { value = std::log(v)/std::log(2.0); }
		logspace( double v ) { value = std::log(v)/std::log(2.0); }
		logspace( long double v ) { value = std::log(v)/std::log(2.0); }
		logspace( const logspace &v ) { value = v.value; }

		// Cast back to regular double.
		inline operator double() const {
			return (double)pow(2.0,value);
		}

		// Cast back to regular float.
		inline operator float() const {
			return (float)pow(2.0,value);
		}

		// Set the log-space value.
		inline logspace &set_base( FPT v ) { value = v; return *this; }

		// Get the log-space value.
		inline FPT as_base() const { return value; }

		// Assignment operators for both ldoubles and from regular double.
		inline logspace &operator=( logspace const &b ) {
			value = b.value;
			return *this;
		}

		inline logspace &operator=( FPT v ) {
			value = std::log(v)/std::log(2.0);
			return *this;
		}

		// Addition uses a lookup table plus a transform.
		inline logspace operator+( logspace b ) const {
			using namespace std;

			if( -value == numeric_limits<FPT>::infinity() ) return b;
			if( -b.value == numeric_limits<FPT>::infinity() ) return *this;

			if( b.value > value ) {
				FPT x = b.value - value;
				if( x > logadd_cutoff ) return b;
				x *= scale;
				b.value += logadd[(unsigned int)x];
				return b;
			}
			// else
			FPT x = value - b.value;
			if( x >= logadd_cutoff ) return *this;
			x *= scale;
			b.value = value + logadd[(unsigned int)x];
			return b;
		}

		inline logspace &operator+=( logspace const &b ) {
			using namespace std;

			if( -value == numeric_limits<FPT>::infinity() ) {
				value = b.value; return *this;
			}
			if( -b.value == numeric_limits<FPT>::infinity() )  return *this;

			if( b.value > value ) {
				value = b.value - value;
				if( value > logadd_cutoff ) { value = b.value; return *this; }
				value *= scale;
				value = logadd[(unsigned int)value];
				value += b.value;
				return *this;
			}
			// else
			FPT x = value - b.value;
			if( x > logadd_cutoff ) return *this;
			x *= scale;
			value += logadd[(unsigned int)x];
			return *this;
		}

		// Multitplication is just addition in log-space.
		inline logspace operator*( logspace v ) const {
			v.value += value;
			return v;
		}

		inline logspace &operator*=( logspace const v ) {
			value += v.value;
			return *this;
		}

		// Division is just subtraction.
		inline logspace operator/( logspace  v ) const {
			v.value = value - v.value;
			return v;
		}

		inline logspace &operator/=( logspace const v ) {
			value -= v.value;
			return *this;
		}

		// Equality and comparison operators do not change.
		inline bool operator==( logspace const v ) const {
			return value == v.value;
		}

		inline bool operator!=( logspace const v ) const {
			return value != v.value;
		}

		inline bool operator<( logspace const v ) const {
			return value < v.value;
		}

		inline bool operator<=( logspace const v ) const {
			return value <= v.value;
		}

		inline bool operator>( logspace const v ) const {
			return value > v.value;
		}

		inline bool operator>=( logspace const v ) const {
			return value >= v.value;
		}

		// Equality and comparison operators to resolve overloading with FPT types.
		inline bool operator==(  FPT v ) const {
			return value == (std::log(v)/std::log(2.0));
		}

		inline bool operator!=(  FPT v ) const {
			return value != (std::log(v)/std::log(2.0));
		}

		inline bool operator<( FPT v ) const {
			return value < (std::log(v)/std::log(2.0));
		}

		inline bool operator<=(  FPT v ) const {
			return value <= (std::log(v)/std::log(2.0));
		}

		inline bool operator>( FPT v ) const {
			return value > (std::log(v)/std::log(2.0));
		}

		inline bool operator>=( FPT v ) const {
			return value >= (std::log(v)/std::log(2.0));
		}

		// Print untransformed value.
		inline friend std::ostream &operator<<( std::ostream &out, logspace const &v ) {
			out << v.value;
			return out;
		}

		//friend logspace pow( logspace v, FPT const x );
	};
	/*
	template< typename FPT >
	FPT logspace<FPT>::logadd[1048576] = {};
	template< typename FPT >
	unsigned int const logspace<FPT>::bits = 20;
	template< typename FPT >
	unsigned int const logspace<FPT>::size = 1048576;
	template< typename FPT >
	FPT const logspace<FPT>::logadd_cutoff = 50.0;
	template< typename FPT >
	FPT const logspace<FPT>::scale = 1048576.0/50.0;
	*/
	template< typename FPT >
	FPT logspace<FPT>::logadd[16777216] = {};
	template< typename FPT >
	unsigned int const logspace<FPT>::bits = 24;
	template< typename FPT >
	unsigned int const logspace<FPT>::size = 16777216;
	template< typename FPT >
	FPT const logspace<FPT>::logadd_cutoff = 50.0;
	template< typename FPT >
	FPT const logspace<FPT>::scale = 16777216.0/50.0;


	template< typename FPT >
	void logspace<FPT>::init() {
		using namespace std;
		FPT x = logadd_cutoff/size;
		for( unsigned int i = 0; i < size; ++i ) {
			FPT b = -(FPT)i*x;
			logadd[i] = log(1.0+pow(2.0,b))/log(2.0);
		}
	}

	// Raising to the power of a FPT is easy, just multiply.
	/*template<typename FPT>
	 inline logspace<FPT> pow( logspace<FPT> v, FPT const x ) {
	 v.value *= x;
	 return v;
	 }*/

};

#endif

