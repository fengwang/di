#ifndef _CRYSTAL_HPP_INCLUDED_FOSIDJ4YOUSHOULDHAVERECEIVEDEOIJSFLKJACOPYOFUTHELCREATEDDIJDOTOWITHDCRYSTALSTRUCT
#define _CRYSTAL_HPP_INCLUDED_FOSIDJ4YOUSHOULDHAVERECEIVEDEOIJSFLKJACOPYOFUTHELCREATEDDIJDOTOWITHDCRYSTALSTRUCT

#include <matrix.hpp>
#include <tri_ary.hpp>

#include <misc/direct_metric_tensor.hpp>
#include <misc/reciprocal_metric_tensor.hpp>
#include <misc/direct_structure_matrix.hpp>
#include <misc/reciprocal_structure_matrix.hpp>
#include <misc/volume.hpp>

namespace feng
{
    
    template<typename T = double>
    struct crystal
    {
        typedef T                       value_type;    
        typedef crystal                 self_type;
        typedef matrix<value_type>      matrix_type;
        typedef tri_ary<value_type>     vector_type;

        //lattice parameters
        value_type a;
        value_type b;
        value_type c;
        value_type alpha;
        value_type beta;
        value_type gamma;

        crystal( const value_type a_, const value_type b_, const value_type c_,
                 const value_type alpha_, const value_type beta_, const value_type gamma_)
            : a(a_), b(b_), c(c_), alpha(alpha_), beta(beta_), gamma(gamma_)
        {}

        crystal( const self_type& ) = default;
        crystal( self_type&& ) = default;
        self_type& operator = ( const self_type& ) = default;
        self_type& operator = ( self_type&& ) = default;

        //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        // vector cross product, dot product, length and included angles
        //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

        value_type
        cartesian_space_cross_product( const vector_type& lhs, const vector_type& rhs ) const 
        {
            return feng::cross_product( lhs, rhs ); 
        }

        value_type
        cartesian_space_dot_product( const vector_type& lhs, const vector_type& rhs ) const 
        {   //implemented in file tri_ary.hpp
            return feng::dot_product(lhs,rhs);        
        }

        value_type
        cartesian_space_length( const vector_type& v ) const 
        {
            return std::sqrt(cartesian_space_dot_product(v,v));
        }

        value_type
        cartesian_space_included_angle( const vector_type& lhs, const vector_type& rhs ) const 
        {
            return std::acos(cartesian_space_dot_product(lhs,rhs)/
                             (cartesian_space_length(lhs)*cartesian_space_length(rhs)));
        }

        //
        value_type
        direct_space_cross_product( const vector_type& lhs, const vector_type& rhs ) const 
        {
            return feng::cross_product( lhs, rhs ) * volume(); 
        }

        value_type 
        direct_space_dot_product( const vector_type& lhs, const vector_type& rhs ) const 
        {   
            return dot_product( lhs, direct_metric_tensor()*rhs );
        }

        value_type
        direct_space_length( const vector_type& v ) const 
        {
            return std::sqrt(direct_space_dot_product(v,v));
        }

        value_type
        direct_space_included_angle( const vector_type& lhs, const vector_type& rhs ) const 
        {
            return std::acos(direct_space_dot_product(lhs,rhs)/
                             (direct_space_length(lhs)*direct_space_length(rhs)));
        }
        
        //
        value_type
        reciprocal_space_cross_product( const vector_type& lhs, const vector_type& rhs ) const 
        {
            return feng::cross_product( lhs, rhs ) / volume(); 
        }

        value_type 
        reciprocal_space_dot_product( const vector_type& lhs, const vector_type& rhs ) const 
        {   
            return dot_product( lhs, reciprocal_metric_tensor()*rhs );
        }

        value_type
        reciprocal_space_length( const vector_type& v ) const 
        {
            return std::sqrt(reciprocal_space_dot_product(v,v));
        }

        value_type
        reciprocal_space_included_angle( const vector_type& lhs, const vector_type& rhs ) const 
        {
            return std::acos(reciprocal_space_dot_product(lhs,rhs)/
                             (reciprocal_space_length(lhs)*reciprocal_space_length(rhs)));
        }

        //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        // cell volume
        //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

        value_type 
        volume() const 
        {   // implemented in file misc/volume.hpp
            return feng::volume(a,b,c,alpha,beta,gamma);
        }

        //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        // direct and reciprocal metric tensors and structure matrice
        //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

        const matrix_type 
        direct_structure_matrix() const 
        {   // implemented in file misc/direct_structure_matrix.hpp
            return std::move(feng::direct_structure_matrix(a,b,c,alpha,beta,gamma));
        }

        const matrix_type
        direct_metric_tensor() const 
        {   // implemented in file misc/direct_metric_tensor.hpp
            return std::move(feng::direct_metric_tensor(a,b,c,alpha,beta,gamma));
        }

        const matrix_type
        reciprocal_structure_matrix() const 
        {   // implemented in file misc/reciprocal_structure_matrix.hpp
            return std::move(feng::reciprocal_structure_matrix(a,b,c,alpha,beta,gamma));
        }

        const matrix_type
        reciprocal_metric_tensor() const 
        {   // implemented in file misc/reciprocal_metric_tensor.hpp
            return std::move(feng::reciprocal_metric_tensor(a,b,c,alpha,beta,gamma));
        }
        
        //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        // vector convertion from onespace to another space
        //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        
        const vector_type
        direct_to_cartesian( const vector_type& d ) const 
        {
            return std::move(direct_structure_matrix() * d);
        }

        const vector_type
        direct_to_reciprocal( const vector_type& d ) const 
        {
            return std::move(d * direct_metric_tensor());
        }
        
        const vector_type
        reciprocal_to_cartesian( const vector_type& r ) const 
        {
            return std::move(reciprocal_structure_matrix() * r);
        }

        const vector_type
        reciprocal_to_direct( const vector_type& r ) const 
        {
            return std::move(r * reciprocal_metric_tensor());
        }

        const vector_type
        cartesian_to_direct( const vector_type& c ) const 
        {
            return std::move(reciprocal_structure_matrix() * c);
        }

        const vector_type
        cartesian_to_reciprocal( const vector_type c ) const 
        {
            return std::move(c * direct_structure_matrix());
        }

    };//crystal

}//namespace feng

#endif//_CRYSTAL_HPP_INCLUDED_FOSIDJ4YOUSHOULDHAVERECEIVEDEOIJSFLKJACOPYOFUTHELCREATEDDIJDOTOWITHDCRYSTALSTRUCT

