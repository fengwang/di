/// @file cartesian_coordinate.h 
///	@brief template class cartesian_coordinate definition and relative operations
/// @author wang feng

#ifndef _COORDINATE_HPP_INCLUDED_SOIJ348YAFOKJ3O8IUWOIUJASFLKJ438UYAFIJFSDOIUJEK
#define _COORDINATE_HPP_INCLUDED_SOIJ348YAFOKJ3O8IUWOIUJASFLKJ438UYAFIJFSDOIUJEK

#include <iostream>
#include <cmath>

namespace feng
{
	//==========================================================================
	//cartesian_coordinate
	//==========================================================================
	/// template calss Cartesian definition.
	/// @param T type of cartesian_coordinate's items
    template < typename T = double >
    struct cartesian_coordinate
    {
			/// @name Standard type definition.

			//@{

			/// The item type, it is T.
            typedef T 						value_type;

			/// The self type, cartesian_coordinate<T>.
            typedef cartesian_coordinate 	self_type;
            
			/// Lvalue of item type. 
			typedef T& 						reference_type;

			/// Const lvalue of item type. 
            //typedef const reference_type 	const_reference_type;
            typedef value_type 				const_reference_type;

			//@}
			
			/// @name Private items 
			/// three dims variables

			//@{

            value_type x_; 			/// < x-axis variable
            value_type y_; 			/// < y-axis variable
            value_type z_; 			/// < z-axis variable

			//@}


			/// @name Constructors.

			//@{

			/// Create from three variables.
            explicit cartesian_coordinate(	const value_type& x = value_type( 0 ),
											const value_type& y = value_type( 0 ),
											const value_type& z = value_type( 0 )
                                         ) : x_( x ), y_( y ), z_( z )
            {}

			/// Create from another cartesian_coordinate.
            cartesian_coordinate( const self_type& other )
                : x_( other.x_ ), y_( other.y_ ), z_( other.z_ )
            {}

			/// Crate from a different type cartesian_coordinate.
            template< typename U >
            cartesian_coordinate( const cartesian_coordinate<U>& other )
                : 	x_( static_cast<T>( other.x_ ) ),
                    y_( static_cast<T>( other.y_ ) ),
                    z_( static_cast<T>( other.z_ ) )
            {}


			//@}


			/// @name Assignments.

			//@{


			/// Assign from another cartesian_coordinate.
            self_type&
            operator = ( const self_type& other )
            {
                x_ = other.x_;
                y_ = other.y_;
                z_ = other.z_;
                return *this;
            }


			/// Assign from a different type cartesian_coordinate.
            template< typename U>
            self_type&
            operator = ( const cartesian_coordinate<U>& other )
            {
                x_ = static_cast<U>( other.x_ );
                y_ = static_cast<U>( other.y_ );
                z_ = static_cast<U>( other.z_ );
                return *this;
            }

			//@}
	



			/// @name Element access.

			//@{

			/// access to x-axis variable, non-const
            reference_type
            x()
            {
                return x_;
            }

			/// access to y-axis variable, non-const
            reference_type
            y()
            {
                return y_;
            }

			/// access to z-axis variable, non-const
            reference_type
            z()
            {
                return z_;
            }

			/// access to x-axis variable, const
			const_reference_type
            x() const
            {
                return x_;
            }

			/// access to y-axis variable, const
            const_reference_type
            y() const
            {
                return y_;
            }

			/// access to z-axis variable, const
            const_reference_type
            z() const
            {
                return z_;
            }

			//@}

			/// Standard output routine for cartesian_coordinate.
			friend
			std::ostream&
			operator << ( std::ostream& lhs, const self_type& rhs )
			{
				lhs << "("
					<< rhs.x_ << ", "
					<< rhs.y_ << ", "
					<< rhs.z_
					<< ")";
				return lhs;
			}

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


			value_type
			square_norm() const
			{
				return x_ * x_ + y_ * y_ + z_ * z_;
			}



	/// Square distance between two cartesian_coordinate positions.

    };//struct cartesian_coordinate



	/// @name  Relative operations.

	//@{

    template< typename T>
    inline
    T
    square_distance( const cartesian_coordinate<T>& lhs, const cartesian_coordinate<T>& rhs )
    {
        const T x_diff = lhs.x() - rhs.x();
        const T y_diff = lhs.y() - rhs.y();
        const T z_diff = lhs.z() - rhs.z();
        const T ans = x_diff * x_diff + y_diff * y_diff + z_diff * z_diff;
        return ans;
    }

	/// Distance between two cartesian_coordinate positions.
    template< typename T>
    inline
    T
    distance( const cartesian_coordinate<T>& lhs, const cartesian_coordinate<T>& rhs )
    {
        const T sd = square_distance( lhs, rhs );
        return std::sqrt( sd );
    }

	/// Equality comparison between two cartesian_coordinate positions.
    template< typename T>
    inline
    bool
    operator == ( const cartesian_coordinate<T>& lhs, const cartesian_coordinate<T>& rhs )
    {
        return 	lhs.x() == rhs.x() &&
                lhs.y() == rhs.y() &&
                lhs.z() == rhs.z();
    }

	/// Unequality comparison between two cartesian_coordinate positions.
    template< typename T>
    inline
    bool
    operator != ( const cartesian_coordinate<T>& lhs, const cartesian_coordinate<T>& rhs )
    {
        return !( lhs == rhs );
    }

	//@}

}//namespace feng

#endif//_COORDINATE_HPP_INCLUDED_SOIJ348YAFOKJ3O8IUWOIUJASFLKJ438UYAFIJFSDOIUJEK

