//
// Adaptor for passing two agents to an observer
//
// by Anton Crombach, A.B.M.Crombach at uu.nl
//

#ifndef _DORSALFIN_DUO_AGENT_
#define _DORSALFIN_DUO_AGENT_

#include "defs.hh"

namespace dorsalfin {

    /// \class DuoAgent
    /// 
    /// Auxilary class for reporting both daughter cells to a logger.
    struct DuoAgent : public Subject, std::pair< Agent*, Agent * > {
        DuoAgent() : Subject(), std::pair< Agent*, Agent* >() {};
        DuoAgent( Agent *a, Agent *b ) 
            : Subject(), std::pair< Agent*, Agent* >( a, b ) {};
        ~DuoAgent() {};
    };
    
    /// For ease of use a make-function
    inline DuoAgent make_duo( Agent *a, Agent *b ) 
    { return DuoAgent( a, b ); }
}
#endif
