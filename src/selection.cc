//
// Implementation of selection schemes.
//
// by Anton Crombach, A.B.M.Crombach at uu.nl
//

#include "selection.hh"

dorsalfin::Agent* 
dorsalfin::RandSelection::select(
        std::vector< Agent* > &va, std::vector< double > &sc ) {
    // note: ignoring scores
    return va[ static_cast< int >( rand_range( va.size() ) ) ];
}

dorsalfin::Agent* 
dorsalfin::ProbalisticSelection::select( 
        std::vector< Agent* > &va, std::vector< double > &sc ) {
    // note: the "no-selection" case may be a 0-pointer in 'va'
    std::partial_sum( sc.begin(), sc.end(), sc.begin() );
    std::vector< double >::iterator aux = 
        std::upper_bound( sc.begin(), sc.end(), rand_range( sc.back() ) );
    int i = std::distance( sc.begin(), aux );
    return va[ i ];
}

