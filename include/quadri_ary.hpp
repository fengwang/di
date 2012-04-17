#ifndef _QUADRI_ARY_HPP_INCLUDED_OSDFIJO4UJSFLKJO4EIFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFDDDD
#define _QUADRI_ARY_HPP_INCLUDED_OSDFIJO4UJSFLKJO4EIFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFDDDD

namespace feng
{
    template<typename T=double>
    struct quadri_ary
    {
        typedef T               value_type;
        typedef quadri_ary      self_type;

        
        value_type u_;
        value_type v_;
        value_type w_;
        value_type t_;

        explicit quadri_ary( const value_type& u = 0, const value_type& v = 0,
                             const value_type& w = 0, const value_type& t = 0 )
            : u_(u), v_(v), w_(w), t_(t) 
        {}

        quadri_ary( const self_type& other ) = default;
        quadri_ary( self_type&& other ) = default;
        self_type& operator = ( const self_type& other ) = default;
        self_type& operator = ( self_type&& other ) = default;

        value_type  u() const { return u_; }
        value_type& u() { return u_; }
        value_type  v() const { return v_; }
        value_type& v() { return v_; }
        value_type  w() const { return w_; }
        value_type& w() { return w_; }
        value_type  t() const { return t_; }
        value_type& t() { return t_; }
    };



}//namespace feng

#endif//_QUADRI_ARY_HPP_INCLUDED_OSDFIJO4UJSFLKJO4EIFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFDDDD

