#ifndef _DIRECT_STRUCTURE_MATRIX_HPP_INCLUDED_FSODIJO4IRU3KLAJFSAJKFLHSDJKDSFLKSFDJALKJFSALFDKJSFDLKJHSAFLKJFLJDDD
#define _DIRECT_STRUCTURE_MATRIX_HPP_INCLUDED_FSODIJO4IRU3KLAJFSAJKFLHSDJKDSFLKSFDJALKJFSALFDKJSFDLKJHSAFLKJFLJDDD

#include <matrix.hpp>
#include <misc/volume.hpp>

#include <cmath>

namespace feng
{
    // Input:
    //          lattice parameters {a, b, c, alpha, beta, gamma}
    // Output:
    //          the direct structure matrix
    template<typename T>
    const matrix<T> 
    direct_structure_matrix( const T a, const T b, const T c, const T alpha, const T beta, const T gamma )
    {
        const T cg = std::cos(gamma);
        const T sg = std::sin(gamma);
        matrix<T> ans(3, 3);
        ans[0][0] = a;   ans[0][1] = b * cg; ans[0][2] = c * std::cos(beta);
        ans[1][0] = T(); ans[1][1] = b * sg; ans[1][2] = - c * (std::cos(beta)*cg-std::cos(alpha)) / sg;
        ans[2][0] = T(); ans[2][1] = T();    ans[2][2] = volume(a,b,c,alpha,beta,gamma)/(a*b*sg);
        return ans;
    }

}//namespace feng

#endif//_DIRECT_STRUCTURE_MATRIX_HPP_INCLUDED_FSODIJO4IRU3KLAJFSAJKFLHSDJKDSFLKSFDJALKJFSALFDKJSFDLKJHSAFLKJFLJDDD

