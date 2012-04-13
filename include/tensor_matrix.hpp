#ifndef _TENSER_MATRIX_INCLUDED_OSFDJHIOI8USLAFDKJ43980UYLSAKFJDHLKSJFADJKDSLFLKSFDSFDLKJSFLDKJSLKJFJFSDJDFDKDFKJKK
#define _TENSER_MATRIX_INCLUDED_OSFDJHIOI8USLAFDKJ43980UYLSAKFJDHLKSJFADJKDSLFLKSFDSFDLKJSFLDKJSLKJFJFSDJDFDKDFKJKK

#include <matrix.hpp>

#include <cmath>

namespace feng
{

template<typename T>
const matrix<T>
make_tensor_matrix( const T a, const T b, const T c, const T alpha, const T beta, const T gama )
{
   matrix<T> g(3, 3); 
   g[0][0] = a*a;
   g[0][1] = a*b * std::cos(gama);
   g[0][2] = a*c * std::cos(beta);
   g[1][0] = g[0][1];
   g[1][1] = b*b;
   g[1][2] = b*c * std::cos(alpha);
   g[2][0] = g[0][2];
   g[2][1] = g[1][2];
   g[2][2] = c*c;

   return g;
}






}

#endif//_TENSER_MATRIX_INCLUDED_OSFDJHIOI8USLAFDKJ43980UYLSAKFJDHLKSJFADJKDSLFLKSFDSFDLKJSFLDKJSLKJFJFSDJDFDKDFKJKK

