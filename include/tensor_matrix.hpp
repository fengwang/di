#ifndef _TENSER_MATRIX_INCLUDED_OSFDJHIOI8USLAFDKJ43980UYLSAKFJDHLKSJFADJKDSLFLKSFDSFDLKJSFLDKJSLKJFJFSDJDFDKDFKJKK
#define _TENSER_MATRIX_INCLUDED_OSFDJHIOI8USLAFDKJ43980UYLSAKFJDHLKSJFADJKDSLFLKSFDSFDLKJSFLDKJSLKJFJFSDJDFDKDFKJKK

#include <matrix.hpp>

#include <cmath>

namespace feng
{

#if 0
    direct metric tensor matrix

            | a^2           ab cos(gamma) ac cos(beta)  |
        g = | ab cos(gamma) b^2           bc cos(alpha) |
            | ac cos(beta)  bc cos(alpha) c^2           |

#endif

template<typename T>
const matrix<T>
tensor_matrix( const T a, const T b, const T c, const T alpha, const T beta, const T gama )
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

}//namespace feng

#endif//_TENSER_MATRIX_INCLUDED_OSFDJHIOI8USLAFDKJ43980UYLSAKFJDHLKSJFADJKDSLFLKSFDSFDLKJSFLDKJSLKJFJFSDJDFDKDFKJKK

