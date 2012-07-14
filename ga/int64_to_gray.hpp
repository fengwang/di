#ifndef _INT64_TO_GRAY_HPP_INCLUDED_FSOIJ3498YUAOFSIUHOUHWAFI8HUEFOHU8ASDVKJBVDIUHAWFGIUYWEFHIUASFIU
#define _INT64_TO_GRAY_HPP_INCLUDED_FSOIJ3498YUAOFSIUHOUHWAFI8HUEFOHU8ASDVKJBVDIUHAWFGIUYWEFHIUASFIU

#include <cstdint>

namespace ga
{
    struct int64_to_gray
    {
        typedef std::uint64_t value_type;

        const value_type operator()( const value_type v ) const 
        {
            return v ^ ( v >> 1 ); 
        }

    };//int64_to_gray
};//namespace ga

#endif//_INT64_TO_GRAY_HPP_INCLUDED_FSOIJ3498YUAOFSIUHOUHWAFI8HUEFOHU8ASDVKJBVDIUHAWFGIUYWEFHIUASFIU

