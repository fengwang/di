#ifndef _DELTA_HPP_INCLUDED_OFDSIJIOURKSLFJ4890ULSAFIJLKFJSDVMNSDFKLIJHFASODIEOIUFSD9084OIJFSDOIU498YFSOIUASFOIUE4
#define _DELTA_HPP_INCLUDED_OFDSIJIOURKSLFJ4890ULSAFIJLKFJSDVMNSDFKLIJHFASODIEOIUFSD9084OIJFSDOIU498YFSOIUASFOIUE4

#include <matrix.hpp>

#include <algorithm>

namespace feng
{
#if 0
                  | 1 0 0 |
    delta_{i,j} = | 0 1 0 |
                  | 0 0 1 |
#endif 
template<typename T>
matrix<T> delta()
{
    matrix<T> ans(3,3);
    std::fill( ans.diag_begin(), ans.diag_end(), T(1) );
    return ans;
}

}//namespace feng

#endif//_DELTA_HPP_INCLUDED_OFDSIJIOURKSLFJ4890ULSAFIJLKFJSDVMNSDFKLIJHFASODIEOIUFSD9084OIJFSDOIU498YFSOIUASFOIUE4

