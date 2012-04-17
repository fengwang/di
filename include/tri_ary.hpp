/// @file tri_ary.hpp
/// @brief Template class tri_ary difinition and relative operations implemented.
/// @author wang feng

#ifndef _TRI_ARY_HPP_INCLUDED_SOIJ349I8AFDSKLJ3OIUAFS98Y34UIOWRE98Y34IOUYHFASDUJ
#define _TRI_ARY_HPP_INCLUDED_SOIJ349I8AFDSKLJ3OIUAFS98Y34UIOWRE98Y34IOUYHFASDUJ

#include <cartesian_coordinate.hpp>
#include <matrix.hpp>

#include <algorithm>
#include <iostream>
#include <cmath>
#include <vector>

namespace feng
{

	/// Template class tri_ary definition.
	/// @param T item type of template class tri_ary 
    template < typename T = double >
    struct tri_ary
    {

			/// @name Public type definition.

			//@{

			/// Item type, T.
            typedef T 						value_type;

			/// Self type, tri_ary<T>.
            typedef tri_ary 				self_type;

			/// Lvalue of items, non-const. 
            typedef value_type&				value_reference;

			// Const lvalue of items.
            typedef const value_reference 	const_value_reference;

			//@}


			/// @name Private items.

			//@{

			/// Component in x direction.
            value_type x_;

			/// Component in y direction.
            value_type y_;

			/// Component in z direction.
            value_type z_;

			//@}

        
			/// @name Constructors.

			//@{

			/// Constructor generated from three variables.
            explicit tri_ary(	const value_type& x = value_type(),
                                const value_type& y = value_type(),
                                const value_type& z = value_type()
                            ) : x_( x ), y_( y ), z_( z )
            {}

			/// Constructor generated from another tri_ary.
            tri_ary( const self_type& other ) :
					x_( other.x_ ),
					y_( other.y_ ),
					z_( other.z_ )
            {}

			/// Constructor generated from a different type tri_ary.
            template< typename U >
            tri_ary( const tri_ary<U>& other ) :
					x_( static_cast<value_type>( other.x_ ) ),
					y_( static_cast<value_type>( other.y_ ) ),
					z_( static_cast<value_type>( other.z_ )
					  )
            {}

			/// Constructor generated from two cartesian_coordinate positions. 
            tri_ary( const cartesian_coordinate<T>& from,
                     const cartesian_coordinate<T>& to )
            {
				x_ = from.x() - to.x();
				y_ = from.y() - to.y();
				z_ = from.z() - to.z();
			}

			//@}


			/// @name Assignments.


			//@{


			/// Assignment from another tri_ary.
            self_type&
            operator = ( const self_type& other )
            {
                x_ = static_cast<value_type>( other.x_ );
                y_ = static_cast<value_type>( other.y_ );
                z_ = static_cast<value_type>( other.z_ );
                return *this;
            }

			/// Assignment from a different type tri_ary.
            template< typename U >
            self_type&
            operator = ( const tri_ary<U>& other )
            {
                x_ = static_cast<value_type>( other.x_ );
                y_ = static_cast<value_type>( other.y_ );
                z_ = static_cast<value_type>( other.z_ );
                return *this;
            }

			//@}

			
			/// @name Unary Operators.

			//@{

			
			/// Times by a same type.
            const self_type&
            operator *=( const value_type& v )
            {
                x_ *= v;
                y_ *= v;
                z_ *= v;
                return *this;
            }

			/// Times by a differenct type.
            template< typename U >
            const self_type&
            operator *=( const U& v )
            {
                x_ *= static_cast<value_type>( v );
                y_ *= static_cast<value_type>( v );
                z_ *= static_cast<value_type>( v );
                return *this;
            }

			/// Divides by a same type.
            const self_type&
            operator /=( const value_type& v )
            {
                x_ /= v;
                y_ /= v;
                z_ /= v;
                return *this;
            }

			/// Divides by a different type.
            template< typename U >
            const self_type&
            operator /=( const U& v )
            {
                x_ /= static_cast<value_type>( v );
                y_ /= static_cast<value_type>( v );
                z_ /= static_cast<value_type>( v );
                return *this;
            }


			/// Plus by another tri_ary.
			/// Implicit conversition allowed.
            const self_type&
            operator += ( const self_type& other )
            {
                x_ += other.x_;
                y_ += other.y_;
                z_ += other.z_;
                return *this;
            }

			/// Minus by another tri_ary.
			/// Implicit conversition allowed.
            const self_type&
            operator -= ( const self_type& other )
            {
                x_ -= other.x_;
                y_ -= other.y_;
                z_ -= other.z_;
                return *this;
            }


			//@}

        


			/// @name Item access.

			//@{

			/// Access to x direction Component, non-const.
            value_reference
            x()
            {
                return x_;
            }

			/// Access to y direction Component, non-const.
            value_reference
            y()
            {
                return y_;
            }

			/// Access to z direction Component, non-const.
            value_reference
            z()
            {
                return z_;
            }

			/// Access to x direction Component, const.
            value_type
            x() const
            {
                return x_;
            }

			/// Access to y direction Component, const.
            value_type
            y() const
            {
                return y_;
            }

			/// Access to z direction Component, const.
            value_type
            z() const
            {
                return z_;
            }


			//@}


			/// @name Functions.


			//@{

			/// Return Euclidean distance.
            value_type
            norm() const
            {
				const value_type m = std::max( 	x_,
												std::max(y_, z_)
											);
				if ( value_type() == m )
					return m;


                const value_type  x = x_ / m;
                const value_type  y = y_ / m;
                const value_type  z = z_ / m;
                const value_type  ans = ::std::sqrt( x*x + y*y + z*z ) * m;

                return ans > value_type() ? ans : -ans;
            }

			/// Return unit tri_ary.
            const self_type
            to_unit() const
            {
                const value_type d = norm();

                if ( value_type() == d )
                {
                    return self_type();
                }

                self_type ans( *this );
                ans /= d;
                return ans;
            }

            // Return a vector with length 3
            const std::vector<value_type>
            to_vector() const
            {
                std::vector<value_type> ans(3);
                ans[0] = x_;
                ans[1] = y_;
                ans[2] = z_;
                return ans;
            }


			//@}
    };


	/// @name binary-operato overloadings

	//@{

    /// Overload operator +.
    template< typename T >
    inline
    const tri_ary<T>
    operator + ( const tri_ary<T>& lhs, const tri_ary<T>& rhs )
    {
        tri_ary<T> ans( lhs );
        ans += rhs;
        return ans;
    }

    /// Overload operator -.
    template< typename T >
    inline
    const tri_ary<T>
    operator - ( const tri_ary<T>& lhs, const tri_ary<T>& rhs )
    {
        tri_ary<T> ans( lhs );
        ans -= rhs;
        return ans;
    }

    /// Overload operator *.
    template< typename T, typename U >
    inline
    const tri_ary<T>
    operator *( const tri_ary<T>& lhs, const U& rhs )
    {
        tri_ary<T> ans( lhs );
        ans *= rhs;
        return ans;
    }

    /// Overload operator /.
    template< typename T, typename U >
    inline
    const tri_ary<T>
    operator / ( const tri_ary<T>& lhs, const U& rhs )
    {
        tri_ary<T> ans( lhs );
        ans /= rhs;
        return ans;
    }

    /// Overload operator <<.
    template<typename T>
    inline
    std::ostream&
    operator << ( std::ostream& lhs, const tri_ary<T>& rhs )
    {
        const char* const tagy = rhs.y() < 0 ? " " : " +";
        const char* const tagz = rhs.z() < 0 ? " " : " +";
        lhs << rhs.x() << "<i>" << tagy << rhs.y() << "<j>" << tagz << rhs.z() << "<k>";
        return lhs;
    }

    /// Overload operator ==.
    template<typename T>
    inline
    bool
    operator == ( const tri_ary<T>& lhs, const tri_ary<T>& rhs )
    {
        return	( lhs.x() == rhs.x() ) &&
                ( lhs.y() == rhs.y() ) &&
                ( lhs.z() == rhs.z() );
    }


	//@}

	/// @name Numerical Methods.

	//@{


	/// Dot product.
    template<typename T>
    inline
    T
    dot_product( const tri_ary<T>& lhs, const tri_ary<T>& rhs )
    {
        const T x = lhs.x() * rhs.x();
        const T y = lhs.y() * rhs.y();
        const T z = lhs.z() * rhs.z();
        const T ans = x + y + z;
        return ans;
    }

	/// Inner Product.
    template<typename T>
    inline
    T
    inner_product( const tri_ary<T>& lhs, const tri_ary<T>& rhs )
    {
        return dot_product( lhs, rhs );
    }

	/// Scalar Product.
    template<typename T>
    inline
    T
    scalar_product( const tri_ary<T>& lhs, const tri_ary<T>& rhs )
    {
        return dot_product( lhs, rhs );
    }

	/// Cross Product.
    template<typename T>
    inline
    const tri_ary<T>
    cross_product( const tri_ary<T>& lhs, const tri_ary<T>& rhs )
    {
        const T x = lhs.y() * rhs.z() - rhs.y() * lhs.z();
        const T y = -lhs.x() * rhs.z() + rhs.x() * lhs.z();
        const T z = lhs.x() * rhs.y() - rhs.x() * lhs.y();
        const tri_ary<T> ans( x, y, z );
        return ans;
    }

	/// Included angle 
    template<typename T>
    inline
    T
    included_angle( const tri_ary<T>& lhs, const tri_ary<T>& rhs )
    {
        const T lhs_norm = lhs.norm();
        const T rhs_norm = rhs.norm();
        const T dot = dot_product( lhs, rhs );

        if ( T() == dot )
        {
            return T();
        }

            const T ans = ::std::acos( dot / ( lhs_norm * rhs_norm ) );

        return ans;
    }
   
    // norm of array
    template<typename T>
    inline 
    T
    norm( const tri_ary<T>& array )
    {
        return array.norm();
    }

    // to unit array
    template<typename T>
    inline 
    const tri_ary<T>  
    to_unit( const tri_ary<T>& array )
    {
        return array.to_unit();
    }

    // to vector
    template<typename T>
    inline
    std::vector<T>
    to_vector(const tri_ary<T>& array)
    {
        return array.to_vector();
    }

    // make tri_ary from two cartesian coordinate objs
    template<typename T>
    inline 
    const tri_ary<T>
    make_tri_ary( const cartesian_coordinate<T>& from, const cartesian_coordinate<T>& to )
    {
        return tri_ary<T>(from, to);
    }

    // 
    template<typename T, std::size_t N, typename A>
    inline 
    const tri_ary<T>
    operator * ( const tri_ary<T>& lhs, const matrix<T,N,A>& rhs )
    {
        assert( 3 == rhs.row() );
        assert( 3 == rhs.col() );
        return tri_ary<T>(  lhs.x() * rhs[0][0] + lhs.y() * rhs[1][0] + lhs.z() * rhs[2][0],
                            lhs.x() * rhs[0][1] + lhs.y() * rhs[1][1] + lhs.z() * rhs[2][1],
                            lhs.x() * rhs[0][2] + lhs.y() * rhs[1][2] + lhs.z() * rhs[2][2] );
    }

    template<typename T, std::size_t N, typename A>
    inline 
    const tri_ary<T>
    operator * ( const matrix<T,N,A>& lhs, const tri_ary<T>& rhs )
    {
        assert( 3 == rhs.row() );
        assert( 3 == rhs.col() );
        return tri_ary<T>(  lhs[0][0] * rhs.x() + lhs[0][1] * rhs.y() + lhs[0][2] * rhs.z(),    
                            lhs[1][0] * rhs.x() + lhs[1][1] * rhs.y() + lhs[1][2] * rhs.z(),    
                            lhs[2][0] * rhs.x() + lhs[2][1] * rhs.y() + lhs[2][2] * rhs.z() );
    }

    template<typename T>
    inline
    T
    operator * ( const tri_ary<T>& lhs, const tri_ary<T>& rhs )
    {
        return dot_product( lhs, rhs );
    }

	//@}

}//namespace feng

#endif//_TRI_ARY_HPP_INCLUDED_SOIJ349I8AFDSKLJ3OIUAFS98Y34UIOWRE98Y34IOUYHFASDUJ

