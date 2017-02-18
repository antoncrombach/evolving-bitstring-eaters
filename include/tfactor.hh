//
// Transcription Factor, a special kind of gene
//
// by Anton Crombach, A.B.M.Crombach@bio.uu.nl
//

#ifndef _DORSALFIN_TFACTOR_H_
#define _DORSALFIN_TFACTOR_H_

#include "defs.hh"
#include "pool.hh"
#include "chromelement.hh"
#include "gene.hh"

namespace dorsalfin {
    
    /// \class TransFactor
    /// \brief Special kind of gene, influences other genes
    class TransFactor : public Gene {
        public:
        /// Ctor
        TransFactor();
        /// Ctor with tag
        TransFactor( int );
        /// Ctor with tag, state and threshold
        TransFactor( int, int, int );
        /// Copy ctor
        explicit TransFactor( const TransFactor & );
        /// Virtual dtor
        virtual ~TransFactor() {}
        // Cloning.
        virtual ChromosomeElement* clone() const;
        /// Return to pool method.
        virtual bool toPool();

        virtual void writeXml( std::ostream & ) const;        
        //
        // the rest is inherited from Gene
        //
    };        
}
#endif
