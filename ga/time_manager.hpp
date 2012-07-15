#ifndef _TIME_MANAGER_HPP_INCLUDED_SDOFJ3498YASFOIUHASFO8UH4O87AGHSFDIUAHSFHTUAWIFOUHASFIUHASFIUHSFI
#define _TIME_MANAGER_HPP_INCLUDED_SDOFJ3498YASFOIUHASFO8UH4O87AGHSFDIUAHSFHTUAWIFOUHASFIUHASFIUHSFI

#include <chrono>
#include <cstdint>

namespace ga
{
    //usage:
    // auto& tm = feng::singleton<time_manager>::instance();
    // tm.set_planning_time( 978987987 ); //within ga manager
    // auto et = tm.elapse_time();
    // auto pt = tm.planning_time();
    // if ( tm.is_timeout() ) ...
    struct time_manager
    {
        std::chrono::time_point<std::chrono::system_clock> start;

        std::uint32_t planning_time_in_second;
        std::uint32_t elapse_time_in_second;

        time_manager( const std::uint32_t planning_time_in_second_ = 1000 ) :
            start( std::chrono::system_clock::now() ),
            planning_time_in_second( planning_time_in_second_ )
        {}

        //set planning time to another value
        void set_planning_time( const std::uint32_t planning_time_in_second_ = 1000 )
        {
            planning_time_in_second = planning_time_in_second_;
        }

        //return time elapsed in second
        std::uint32_t elapse_time()
        {
            auto const current_time = std::chrono::system_clock::now();
            elapse_time_in_second = std::chrono::duration_cast<std::chrono::seconds> ( current_time - start ).count();
            return elapse_time_in_second;
        }

        std::uint32_t planning_time() const
        {
            return planning_time_in_second;
        }

        bool is_timeout()
        {
            return elapse_time() > planning_time_in_second;
        }

    };//time_manager

}//namespace ga

#endif//_TIME_MANAGER_HPP_INCLUDED_SDOFJ3498YASFOIUHASFO8UH4O87AGHSFDIUAHSFHTUAWIFOUHASFIUHASFIUHSFI

