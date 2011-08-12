#ifndef __BTL_DNA_IUPAC_GAP__
#define __BTL_DNA_IUPAC_GAP__

#include <inttypes.h>

namespace btl {
   class dna_iupac_gap {
      public:
         typedef uint8_t symbol;
         static symbol const alphabet[];
         static symbol const alphabet_decode[];
         static int const symbol_ambiguity[];
         typedef symbol const * const_iterator;

         enum symbols {GAP=0, G, C, S, T, K, Y, B, A, R, M, V, W, D, H, N, end_token};

         inline static int size() { 
            return end_token;
         }

         inline static const_iterator begin() {
            return &alphabet[GAP];
         }
         
         inline static const_iterator end() {
            return &alphabet[end_token];
         }

         inline static symbol consensus( symbol c1, symbol c2 ) {
            return ((c1 & c2) == 0) ? (c1 | c2) : (c1 & c2);
         }

         inline static symbol decode( symbol const b ){
            return alphabet_decode[(int)b];
         }

         inline static symbol encode(symbol const c ) {
            switch( c ) {
               case '-':
                  return GAP;

               case 'g':
               case 'G':
                  return G;

               case 'c':
               case 'C':
                  return C;

               case 's':
               case 'S':
                  return S;

               case 't':
               case 'T':
                  return T;

               case 'k':
               case 'K':
                  return K;

               case 'y':
               case 'Y':
                  return Y;

               case 'b':
               case 'B':
                  return B;

               case 'a':
               case 'A':
                  return A;

               case 'r':
               case 'R':
                  return R;

               case 'm':
               case 'M':
                  return M;

               case 'v':
               case 'V':
                  return V;

               case 'w':
               case 'W':
                  return W;

               case 'd':
               case 'D':
                  return D;

               case 'h':
               case 'H':
                  return H;

               case 'n':
               case 'N':
               case 'x':
               case 'X':
                  return N;

               default:
                  // handle error?
                  return end_token;
            }
         }

         inline static int ambiguity( symbol const b ) {
            return symbol_ambiguity[(int)b];
         }

   };
}

#endif
