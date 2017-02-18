//
// Implementation of events.
//
// by Anton Crombach, A.B.M.Crombach at uu.nl
//

#include "event.hh"
#include "chromosome.hh"

void
dorsalfin::ZeroRatesEvent::execute( Population &pop ) {
    Chromosome::zeroRates();
}
