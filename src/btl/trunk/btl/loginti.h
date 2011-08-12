/** Log space fixed point number.
 *
 *  This class uses fake fixed point numbers in logspace. Normally one would designate a
 *  number of bits to store the fractional part of a fixed point number. Here instead,
 *  we treat the number as base 10 and just assume that the last d digits are the fractional
 *  point. For example, if we use int32 as the base type, the maximum value is
 *  2,147,483,648. If we assume 4 digits of precision, the maximum value becomes
 *  214,748.3648
 */
#ifndef BTL_LOGINTI_H__
#define BTL_LOGINTI_H__

#include <cmath>
#include <limits>
#include <iostream>

namespace btl {

	// static lookup tables
	template< typename BASET, int D >
	class loginti {
		struct table_entry {
			int y;
			int slope;
		};

	public:
		static table_entry *logadd;
		static BASET const logadd_cutoff;
		static unsigned int const table_size;
		static double const increment;
		static double const offset;
		static int const scaleI;
		static double const scaleD;
		static float const scaleF;


		// compute tables
		static void init();

	private:
		BASET value;

	public:

		// Constructors for various types.
		loginti() {}
		loginti( float v ) { value = scaleF*std::log(v)/std::log(2.0); }
		loginti( double v ) { value = scaleD*std::log(v)/std::log(2.0); }
		loginti( long double v ) { value = scaleD*std::log(v)/std::log(2.0); }
		loginti( const loginti &v ) { value = v.value; }

		// Cast back to regular double.
		inline operator double() const {
			return pow(2.0,(double)value/scaleD);
		}

		// Cast back to regular float.
		inline operator float() const {
			return (float)pow(2.0,(float)value/scaleF);
		}

		// Set the log-space value.
		inline loginti &set_base( double v ) { value = scaleD*v; return *this; }

		// Get the log-space value.
		inline double as_base() const { return value/scaleD; }

		// Assignment operators for both ldoubles and from regular double.
		inline loginti &operator=( loginti const &b ) {
			value = b.value;
			return *this;
		}

		inline loginti &operator=( double v ) {
			value = scaleD*std::log(v)/std::log(2.0);
			return *this;
		}

		// Addition uses a lookup table plus a tricky transform.
		inline loginti operator+( loginti b ) const {
			using namespace std;

			if( -value == numeric_limits<BASET>::infinity() ) { return b; cerr << "FAIL" << endl; }
			if( -b.value == numeric_limits<BASET>::infinity() ) { return *this; cerr << "FAIL" << endl; }

			if( b.value > value ) {
				BASET x = b.value - value;
				int index = x >> 8;
				if( x > logadd_cutoff ) return b;
				table_entry entry = logadd[index];
				b.value += entry.y + entry.slope*x/scaleI;
				return b;
			}
			// else
			BASET x = value - b.value;
			int index = x >> 8;
			if( x > logadd_cutoff ) return *this;
			table_entry entry = logadd[index];
			b.value = value + entry.y + entry.slope*x/scaleI;
			return b;
		}

		inline loginti &operator+=( loginti const &b ) {
			using namespace std;
			exit(1);

			if( -value == numeric_limits<BASET>::infinity() ) {
				value = b.value; return *this;
			}
			if( -b.value == numeric_limits<BASET>::infinity() )  return *this;

			if( b.value > value ) {
				value = b.value - value;
				if( value > logadd_cutoff ) { value = b.value; return *this; }
				value = logadd[value/(int)increment].y + logadd[value/(int)increment].slope*value/(int)scaleI;
				value += b.value;
				return *this;
			}
			// else
			BASET x = value - b.value;
			if( x > logadd_cutoff ) return *this;
			value += logadd[x/(int)increment].y + logadd[x/(int)increment].slope*x/scaleI;
			return *this;
		}

		// Multiplication is just addition in log-space.
		inline loginti operator*( loginti v ) const {
			v.value += value;
			return v;
		}

		inline loginti &operator*=( loginti const v ) {
			value += v.value;
			return *this;
		}

		// Division is just subtraction.
		inline loginti operator/( loginti  v ) const {
			v.value = value - v.value;
			return v;
		}

		inline loginti &operator/=( loginti const v ) {
			value -= v.value;
			return *this;
		}

		// Equality and comparison operators do not change.
		inline bool operator==( loginti const v ) const {
			return value == v.value;
		}

		inline bool operator!=( loginti const v ) const {
			return value != v.value;
		}

		inline bool operator<( loginti const v ) const {
			return value < v.value;
		}

		inline bool operator<=( loginti const v ) const {
			return value <= v.value;
		}

		inline bool operator>( loginti const v ) const {
			return value > v.value;
		}

		inline bool operator>=( loginti const v ) const {
			return value >= v.value;
		}

		// Equality and comparison operators to resolve overloading with FPT types.
		inline bool operator==( double v ) const {
			return value == (scaleD*std::log(v)/std::log(2.0));
		}

		inline bool operator!=( double v ) const {
			return value != (scaleD*std::log(v)/std::log(2.0));
		}

		inline bool operator<( double v ) const {
			return value < (scaleD*std::log(v)/std::log(2.0));
		}

		inline bool operator<=( double v ) const {
			return value <= (scaleD*std::log(v)/std::log(2.0));
		}

		inline bool operator>( double v ) const {
			return value > (scaleD*std::log(v)/std::log(2.0));
		}

		inline bool operator>=( double v ) const {
			return value >= (scaleD*std::log(v)/std::log(2.0));
		}

		// Print untransformed value.
		inline friend std::ostream &operator<<( std::ostream &out, loginti const &v ) {
			out << (double)v.value/scaleD;
			return out;
		}


	};

	template< typename BASET, int D >
	double const loginti<BASET,D>::scaleD = std::pow(10.0,(double)D);

	template< typename BASET, int D >
	float const loginti<BASET,D>::scaleF = (float)scaleD;

	template< typename BASET, int D >
	int const loginti<BASET,D>::scaleI = (int)scaleD;


	template< typename BASET, int D >
	BASET const loginti<BASET,D>::logadd_cutoff = (BASET) (-scaleD*std::log(std::pow(2.0,1.0/scaleD) - 1.0)/std::log(2.0));

	template< typename BASET, int D >
	double const loginti<BASET,D>::increment = 256;

	template< typename BASET, int D >
	unsigned int const loginti<BASET,D>::table_size = logadd_cutoff/increment + 1;

	template< typename BASET, int D >
	typename loginti<BASET,D>::table_entry* loginti<BASET,D>::logadd = 0;

	template< typename BASET, int D >
	void loginti<BASET,D>::init() {
		using namespace std;
		logadd = new table_entry[table_size];
		cerr << table_size << '\t' << increment << endl;
		double increment = (double)logadd_cutoff/(double)table_size;
		double j = 0.0, cur, next;
		for( unsigned int i = 0; i <= table_size; ++i, j += increment ) {
			double b = -j/scaleD, c = -(j+increment)/scaleD;
			cur = log(1.0+pow(2.0,b))/log(2.0);
			next = log(1.0+pow(2.0,b))/log(2.0);
			logadd[i].y = scaleD*(cur - ((next-cur)/increment)*j);
			logadd[i].slope = scaleD*((next-cur)/increment);
		}
	}


};

#endif

