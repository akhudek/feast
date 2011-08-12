#ifndef IO_H__
#define IO_H__
#include <iostream>
#include <string>
#include <stdexcept>

// reads a type from a file, skipping whitespace and comments
template<typename MyType>
MyType read( std::istream &in )  {
	while( 1 ) {
      int ch = in.peek();
      if( !in.good() ) throw std::runtime_error("read: failed in file access");
		
      // is a comment, so skip to next line
      if( ch == '#' ) {
         while( in.get() != '\n' ) if( !in.good() ) throw std::runtime_error("read: failed in comment");
         continue;
      }
		
      // skip whitespace
      if( std::isspace(ch) ) { in.get(); continue; }
		
      // assume it's a MyType next
      break;
   }
	
   MyType value;
   in >> value;
   if( !in.good() ) throw std::runtime_error("read: failed reading");
   return value;
}

#endif
