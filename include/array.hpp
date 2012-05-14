#ifndef _ARRAY_HPP_INCLUDED_SDFIOJ9438USLKFDJDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDSFD8E3I
#define _ARRAY_HPP_INCLUDED_SDFIOJ9438USLKFDJDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDSFD8E3I

#include <array>
#include <algorithm>
#include <iosfwd>
#include <cassert>

#include <matrix.hpp>

#include <iostream>

namespace feng
{

    template<typename T>
    struct array : std::array<T,3>
    {
        typedef T               value_type;
        typedef array           self_type;

        value_type x() const { return (*this)[0]; }
        value_type y() const { return (*this)[1]; }
        value_type z() const { return (*this)[2]; }

        value_type& x() { return (*this)[0]; }
        value_type& y() { return (*this)[1]; }
        value_type& z() { return (*this)[2]; }

        explicit array ( const value_type a = 0, const value_type b = 0, const value_type c = 0 )
        {
            (*this)[0] = a; (*this)[1] = b; (*this)[2] = c;
        }
        
        value_type norm() const 
        {
            const value_type m = std::max(std::abs((*this)[0]), 
                                          std::max(std::abs((*this)[1]), std::abs((*this)[2])));
            if ( value_type() == m ) return value_type();
            const value_type x = (*this)[0]/m;
            const value_type y = (*this)[1]/m;
            const value_type z = (*this)[2]/m;
            return m * std::sqrt(x*x +y*y +z*z);
        }

        const self_type operator - () const 
        {
            return self_type(-(*this)[0], -(*this)[1], -(*this)[2]);
        }

        friend std::ostream& operator << ( std::ostream& os, const self_type& self )
        {
            return os <<          self[0] << " "
                               << self[1] << " "
                               << self[2];
            //return os << " ("  << self[0] << ", "
            //                   << self[1] << ", "
            //                   << self[2] << ") " ;
        }

        template<std::size_t D, typename A>
        friend const self_type operator * ( const matrix<T,D,A>& lhs, const self_type& rhs )
        {
            assert( 3 == lhs.row() );
            assert( 3 == lhs.col() );
            return self_type(   std::inner_product( lhs.row_begin(0), lhs.row_end(0), rhs.begin(), value_type() ),
                                std::inner_product( lhs.row_begin(1), lhs.row_end(1), rhs.begin(), value_type() ),
                                std::inner_product( lhs.row_begin(2), lhs.row_end(2), rhs.begin(), value_type() ));
        }

        template<std::size_t D, typename A>
        friend const self_type operator * ( const self_type& lhs, const matrix<T,D,A>& rhs )
        {
            assert( 3 == rhs.row() ); 
            assert( 3 == rhs.col() ); 
            return self_type (  std::inner_product( lhs.begin(), lhs.end(), rhs.col_begin(0), value_type() ),
                                std::inner_product( lhs.begin(), lhs.end(), rhs.col_begin(1), value_type() ),
                                std::inner_product( lhs.begin(), lhs.end(), rhs.col_begin(2), value_type() ));
        }
        
        friend const self_type operator + ( const self_type& lhs, const self_type rhs )
        {
            return self_type(lhs[0]+rhs[0], lhs[1]+rhs[1], lhs[2]+rhs[2]);
        }
        
        friend const self_type operator + ( const self_type& lhs, const value_type rhs )
        {
            return self_type(lhs[0]+rhs, lhs[1]+rhs, lhs[2]+rhs);
        }
        
        friend const self_type operator + ( const value_type lhs, const self_type& rhs )
        {
            return rhs + lhs;
        }
        
        friend const self_type operator - ( const self_type& lhs, const value_type rhs )
        {
            return self_type(lhs[0]-rhs, lhs[1]-rhs, lhs[2]-rhs);
        }
        
        friend const self_type operator - ( const self_type& lhs, const self_type rhs )
        {
            return self_type(lhs[0]-rhs[0], lhs[1]-rhs[1], lhs[2]-rhs[2]);
        }

        friend const self_type operator * ( const self_type& lhs, const value_type rhs )
        {
            return self_type(lhs[0]*rhs, lhs[1]*rhs, lhs[2]*rhs);
        }
        
        friend const self_type operator * ( const value_type lhs, const self_type& rhs )
        {
            return rhs * lhs;
        }
        
        friend const self_type operator / ( const self_type& lhs, const value_type rhs )
        {
            return self_type(lhs[0]/rhs, lhs[1]/rhs, lhs[2]/rhs);
        }
        
        friend bool operator == ( const self_type& lhs, const self_type& rhs )
        {
            return  lhs[0] == rhs[0] &&
                    lhs[1] == rhs[1] &&
                    lhs[2] == rhs[2];
        }
        
        friend bool operator != ( const self_type& lhs, const self_type& rhs )
        {
            return  lhs[0] != rhs[0] ||
                    lhs[1] != rhs[1] ||
                    lhs[2] != rhs[2];
        }
        
        friend bool operator < ( const self_type& lhs, const self_type& rhs )
        {
            if ( lhs[0] < rhs[0] ) return true;
            if ( lhs[0] > rhs[0] ) return false;
            if ( lhs[1] < rhs[1] ) return true;
            if ( lhs[1] > rhs[1] ) return false;
            if ( lhs[2] < rhs[2] ) return true;
            return false;
        }
        
        friend bool operator <= ( const self_type& lhs, const self_type& rhs )
        {
            return (lhs < rhs) || (lhs == rhs);
        }
        
        friend bool operator > ( const self_type& lhs, const self_type& rhs )
        {
            if ( lhs[0] > rhs[0] ) return true;
            if ( lhs[0] < rhs[0] ) return false;
            if ( lhs[1] > rhs[1] ) return true;
            if ( lhs[1] < rhs[1] ) return false;
            if ( lhs[2] > rhs[2] ) return true;
            return false;
        }
        
        friend bool operator >= ( const self_type& lhs, const self_type& rhs )
        {
            return (lhs > rhs) || (lhs == rhs);
        }

    };//struct array

    template<typename T>
    T inner_product( const array<T>& lhs, const array<T>& rhs )
    {
        return std::inner_product(lhs.begin(), lhs.end(), rhs.begin(), T(0));
    }

    template<typename T>
    const array<T> cross_product( const array<T>& lhs, const array<T>& rhs )
    {
        return array<T>( lhs[1]*rhs[2] - rhs[1]*lhs[2],
                         lhs[2]*rhs[0] - rhs[2]*lhs[0],
                         lhs[0]*rhs[1] - rhs[0]*lhs[1] );
    }

    template<typename T>
    T included_angle( const array<T>& lhs, const array<T>& rhs )
    {
        const T ln = lhs.norm();
        const T rn = rhs.norm();
        const T lr = inner_product(lhs, rhs);
        if ( T() == lr ) return T();
        return std::acos(lr/(ln*rn));
    }

};//namespace feng

#endif//_ARRAY_HPP_INCLUDED_SDFIOJ9438USLKFDJDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDSFD8E3I

