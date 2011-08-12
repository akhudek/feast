/*
 *  gff_feature.h
 *  btl
 *
 *  Created by Alexander K. Hudek on 2008-09-17.
 *  Copyright 2008 University of Waterloo. All rights reserved.
 *
 */
#ifndef BTL_GFF_FEATURE_H__
#define BTL_GFF_FEATURE_H__
#include <btl/annotation.h>
#include <boost/lexical_cast.hpp>
#include <stdexcept>
#include <string>
#include <map>

namespace btl {

	class gff_feature : public annotation {
	public:
		enum strand_type { NOSTRAND = 0, POSITIVE, NEGATIVE };
		enum frame_type { F0 = 0, F1 = 1, F2 = 2, FNONE };

		std::string name;

		std::string source;

		std::string feature;

		double score;

		strand_type strand;

		frame_type frame;

		std::map< std::string, std::string > attributes;

		// default constructor
		gff_feature() : annotation() {}

		// read a feature from a gff line
		inline void assign_from_string( std::string const &in ) {
			using namespace std;
			string::size_type prev = 0, pos = 0;

			// get name
			pos = in.find( '\t', prev );
			if( pos == string::npos ) throw runtime_error("gff_feature: missing name");
			name = in.substr( prev, pos - prev );
			prev = pos + 1;

			// get source
			pos = in.find( '\t', prev );
			if( pos == string::npos ) throw runtime_error("gff_feature: missing source");
			source = in.substr( prev, pos - prev );
			prev = pos + 1;

			// get feature
			pos = in.find( '\t', prev );
			if( pos == string::npos ) throw runtime_error("gff_feature: missing feature");
			feature = in.substr( prev, pos - prev );
			prev = pos + 1;

			// get start
			pos = in.find( '\t', prev );
			if( pos == string::npos ) throw runtime_error("gff_feature: missing start");
			start = boost::lexical_cast<int>( in.substr( prev, pos - prev ) );
			prev = pos + 1;

			// get end
			pos = in.find( '\t', prev );
			if( pos == string::npos ) throw runtime_error("gff_feature: missing end");
			end = boost::lexical_cast<int>( in.substr( prev, pos - prev ) );
			prev = pos + 1;

			// get score
			pos = in.find( '\t', prev );
			score = boost::lexical_cast<double>( in.substr( prev, pos - prev ) );
			prev = pos + 1;

			// get strand
			pos = in.find( '\t', prev );
			if( pos == string::npos ) throw runtime_error("gff_feature: missing strand");
			std::string strand_str = in.substr( prev, pos - prev );
			prev = pos + 1;

			if( strand_str.size() != 1 ) throw runtime_error("gff_feature: invalid strand feild");
			switch( strand_str[0] ) {
				case '+':
					strand = POSITIVE;
					break;
				case '-':
					strand = NEGATIVE;
					break;
				case '.':
					strand = NOSTRAND;
					break;
				default:
					throw runtime_error("gff_feature: invalid strand feild");
			}

			// get frame
			pos = in.find( '\t', prev );
			if( pos == string::npos ) throw runtime_error("gff_feature: missing frame");
			std::string frame_str = in.substr( prev, pos - prev );
			prev = pos + 1;

			if( frame_str.size() != 1 ) throw runtime_error("gff_feature: invalid frame field");
			switch( frame_str[0] ) {
				case '0':
					frame = F0;
					break;
				case '1':
					frame = F1;
					break;
				case '2':
					frame = F2;
					break;
				case '.':
					frame = FNONE;
					break;
				default:
					throw runtime_error("gff_feature: invalid frame field");
			}

			// TODO read attributes
		}

		// Read a feature from a gff line. We skip blank lines and comments.
		inline void assign_from_stream( std::istream &in ) {
			std::string instr;
			std::getline(in,instr);
			assign_from_string(instr);
		}

		inline gff_feature &operator=( gff_feature const &o ) {
			name = o.name;
			source = o.source;
			feature = o.feature;
			start = o.start;
			end = o.end;
			type = o.type;
			score = o.score;
			strand = o.strand;
			frame = o.frame;
			attributes = o.attributes;
			return *this;
		}

	};
}

#endif
