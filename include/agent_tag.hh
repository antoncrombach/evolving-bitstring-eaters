//
// Agent identification tag
//
// by Anton Crombach, A.B.M.Crombach at uu.nl
//

#ifndef _DORSALFIN_AGENT_TAG_
#define _DORSALFIN_AGENT_TAG_

#include "defs.hh"

namespace dorsalfin {
    
    /// \class AgentTag
    /// \brief Unique identification for agents
    /// 
    /// Ancestor tracing of agents is done with a small \c class named
    /// \c AgentTag.
    class AgentTag {
        public:
        /// Constructor for convenience
        AgentTag() : time( -1 ), loc( std::make_pair( -1, -1 ) ), i( 0 ) {}
        /// Constructor
        AgentTag( long t, location l, int i ) 
            : time( t ), loc( l ), i( i ) {}
        /// Destructor
        ~AgentTag() {}

        /// Reset all fields to null values
        void reset();
        /// Test for null values
        bool isNull() const;
        
        /// As a string
        std::string str() const;

        public:        
        /// Time of birth
        long time;
        /// Birth location
        location loc;
        /// Index for those agents that are born in the same timestep as their
        /// parent
        int i;
    };
    
    inline bool operator==( const AgentTag &t1, const AgentTag &t2 ) { 
        return t1.time == t2.time and t1.loc == t2.loc and t1.i == t2.i; 
    }
    
    inline bool operator<( const AgentTag &t1, const AgentTag &t2 ) { 
        return t1.time < t2.time or 
            ( t1.time == t2.time and t1.loc < t2.loc ) or    
            ( t1.time == t2.time and t1.loc == t2.loc and t1.i < t2.i );
    }
    
    inline std::ostream & operator<<( std::ostream &os, const AgentTag &t ) {
        os <<  t.time << "\t" << t.loc.first << "\t" << t.loc.second 
           << "\t" << t.i;
        return os;
    }
    
    inline void AgentTag::reset() {
        time = -1;
        loc = std::make_pair( -1, -1 );
        i = 0;
    }
    
    inline bool AgentTag::isNull() const {
        return time == -1 and loc == std::make_pair( -1, -1 ) and i == 0;
    }
    
    inline std::string AgentTag::str() const
    { return boost::lexical_cast< std::string >( time ) + "-" +
        boost::lexical_cast< std::string >( loc.first ) + "-" +
        boost::lexical_cast< std::string >( loc.second ) + "-" +
        boost::lexical_cast< std::string >( i ); }
}
#endif
