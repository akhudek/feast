/** Log space fixed point number.
 *
 *  This class uses fake fixed point numbers in logspace. Normally one would designate a
 *  number of bits to store the fractional part of a fixed point number. Here instead,
 *  we treat the number as base 10 and just assume that the last d digits are the fractional
 *  part. For example, if we use int32 as the base type, the maximum value is
 *  2,147,483,648. If we assume 4 digits of precision, the maximum value becomes
 *  214,748.3648
 */
#ifndef BTL_LOGINT_H__
#define BTL_LOGINT_H__

#include <cmath>
#include <limits>
#include <iostream>

namespace btl {

	// static lookup tables
	template< typename BASET, int D >
	class logint {

	public:
		static BASET *logadd;
		static BASET const logadd_cutoff;
		static double const scaleD;
		static float const scaleF;

		// compute tables
		static void init();

	private:
		BASET value;

	public:

		// Constructors for various types.
		logint() {}
		logint( float v ) { value = scaleF*std::log(v)/std::log(2.0); }
		logint( double v ) { value = scaleD*std::log(v)/std::log(2.0); }
		logint( long double v ) { value = scaleD*std::log(v)/std::log(2.0); }
		logint( const logint &v ) { value = v.value; }

		// Cast back to regular double.
		inline operator double() const {
			return pow(2.0,(double)value/scaleD);
		}

		// Cast back to regular float.
		inline operator float() const {
			return (float)pow(2.0,(float)value/scaleF);
		}

		// Set the log-space value.
		inline logint &set_base( double v ) { value = scaleD*v; return *this; }

		// Get the log-space value.
		inline double as_base() const { return (double)value/scaleD; }

		// Assignment operators for both ldoubles and from regular double.
		inline logint &operator=( logint const &b ) {
			value = b.value;
			return *this;
		}

		inline logint &operator=( double v ) {
			value = scaleD*std::log(v)/std::log(2.0);
			return *this;
		}

		// Addition uses a lookup table plus a tricky transform.
		inline logint operator+( logint b ) const {
			using namespace std;

			if( value == numeric_limits<BASET>::min() ) return b;
			if( b.value == numeric_limits<BASET>::min() ) return *this;

			if( b.value > value ) {
				BASET x = b.value - value;
				if( x > logadd_cutoff ) return b;
				b.value += logadd[x];
				return b;
			}
			// else
			BASET x = value - b.value;
			if( x > logadd_cutoff ) return *this;
			b.value = value + logadd[x];
			return b;
		}

		inline logint &operator+=( logint const &b ) {
			using namespace std;

			if( value == numeric_limits<BASET>::min() ) {
				value = b.value; return *this;
			}
			if( b.value == numeric_limits<BASET>::min() )  return *this;

			if( b.value > value ) {
				BASET x = b.value - value;
				if( x > logadd_cutoff ) { value = b.value; return *this; }
				value = b.value + logadd[x];
				return *this;
			}
			// else
			BASET x = value - b.value;
			if( x > logadd_cutoff ) return *this;
			value += logadd[x];
			return *this;
		}

		// Multiplication is just addition in log-space.
		inline logint operator*( logint v ) const {
			using namespace std;
			if( value == numeric_limits<BASET>::min() || v.value == numeric_limits<BASET>::min() )
				v.value = numeric_limits<BASET>::min();
			else v.value += value;
			return v;
		}

		inline logint &operator*=( logint const v ) {
			using namespace std;
			if( value == numeric_limits<BASET>::min() || v.value == numeric_limits<BASET>::min() )
				value = numeric_limits<BASET>::min();
			else value += v.value;
			return *this;
		}

		// Division is just subtraction.
		inline logint operator/( logint  v ) const {
			using namespace std;
			if( value == numeric_limits<BASET>::min() ) v.value = numeric_limits<BASET>::min();
			else v.value = value - v.value;
			return v;
		}

		inline logint &operator/=( logint const v ) {
			using namespace std;
			if( value != numeric_limits<BASET>::min() ) value -= v.value;
			return *this;
		}

		// Equality and comparison operators do not change.
		inline bool operator==( logint const v ) const {
			return value == v.value;
		}

		inline bool operator!=( logint const v ) const {
			return value != v.value;
		}

		inline bool operator<( logint const v ) const {
			return value < v.value;
		}

		inline bool operator<=( logint const v ) const {
			return value <= v.value;
		}

		inline bool operator>( logint const v ) const {
			return value > v.value;
		}

		inline bool operator>=( logint const v ) const {
			return value >= v.value;
		}

		// Equality and comparison operators to resolve overloading with FPT types.
		inline bool operator==( double v ) const {
			return value == (BASET)(scaleD*std::log(v)/std::log(2.0));
		}

		inline bool operator!=( double v ) const {
			return value != (BASET)(scaleD*std::log(v)/std::log(2.0));
		}

		inline bool operator<( double v ) const {
			return value < (BASET)(scaleD*std::log(v)/std::log(2.0));
		}

		inline bool operator<=( double v ) const {
			return value <= (BASET)(scaleD*std::log(v)/std::log(2.0));
		}

		inline bool operator>( double v ) const {
			return value > (BASET)(scaleD*std::log(v)/std::log(2.0));
		}

		inline bool operator>=( double v ) const {
			return value >= (BASET)(scaleD*std::log(v)/std::log(2.0));
		}

		// Print untransformed value.
		inline friend std::ostream &operator<<( std::ostream &out, logint const &v ) {
			out << (double)v.value/scaleD;
			return out;
		}


	};

	template< typename BASET, int D >
	double const logint<BASET,D>::scaleD = std::pow(10.0,(double)D);

	template< typename BASET, int D >
	float const logint<BASET,D>::scaleF = (float)scaleD;

	template< typename BASET, int D >
	BASET const logint<BASET,D>::logadd_cutoff = (BASET) (-scaleD*std::log(std::pow(2.0,1.0/scaleD) - 1.0)/std::log(2.0));

	template< typename BASET, int D >
	BASET* logint<BASET,D>::logadd = 0;

	/** This initializes lookup tables required by the logint class.
	 */
	template< typename BASET, int D >
	void logint<BASET,D>::init() {
		using namespace std;
		logadd = new BASET[logadd_cutoff+1];
		for( unsigned int i = 0; i <= logadd_cutoff; ++i ) {
			double b = -(double)i/scaleD;
			logadd[i] = scaleD*log(1.0+pow(2.0,b))/log(2.0);
		}
	}


};

#endif

