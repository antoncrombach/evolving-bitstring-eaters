//
// Implementation of abstract agent.
//
// by Anton Crombach, A.B.M.Crombach at uu.nl
//

#include "agent.hh"

dorsalfin::Agent::Agent() 
: dying_( false ), just_born_( true ), me_(), ancestor_(), type_( -1 ) {}

dorsalfin::Agent::Agent( AgentTag t ) 
: dying_( false ), just_born_( true ), me_( t ), ancestor_(), type_( -1 ) {}

dorsalfin::Agent::Agent( int tt ) 
: dying_( false ), just_born_( true ), me_(), ancestor_(), type_( tt ) {}

dorsalfin::Agent::Agent( const Agent &ag ) {
    copy( ag );
}

void
dorsalfin::Agent::copy( const Agent &ag ) {
    dying_ = ag.dying_;
    me_ = ag.me_;
    ancestor_ = ag.ancestor_;
    type_ = ag.type_;
    just_born_ = ag.just_born_;
}

