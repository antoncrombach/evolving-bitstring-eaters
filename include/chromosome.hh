//
// Chromosome of an agent, part of the genome.
//
// by Anton Crombach, A.B.M.Crombach at uu.nl
//

#ifndef _DORSALFIN_CHROMOSOME_H_
#define _DORSALFIN_CHROMOSOME_H_

#include "defs.hh"
#include "chromelement.hh"
#include "genome.hh"
#include "gene.hh"
#include "tfactor.hh"
#include "interaction.hh"
#include "repeat.hh"
#include "retroposon.hh"
#include "distribution.hh"

namespace dorsalfin {
   
    /// \class Chromosome
    /// \brief A sequence of genes, repeats and retroposons.
    ///
    /// The heart of the genome model. Basically, it is a sequence of 
    /// genes, bsites, long terminal repeats and retroposons in any
    /// order. 
    ///
    /// From this sequence, a transcription network can be built.
    class Chromosome : public CachedElement, public XmlWriter {
        public:
        // We had (6) events, now it's 11, though DSB, CP_RP, RM_RP and RM_LTR
        // are not used for the time being.
        enum mut_event { DSB, CP_G, RM_G, THR, TAG, CP_RP, RM_RP, RM_LTR,  
            W_IA, T_IA, CP_IA, RM_IA, NW_IA, ST };
        
        public:
        /// Chromosome iterator
        typedef std::list< ChromosomeElement* >::iterator ce_iter;
        /// Chromosome const iterator
        typedef std::list< ChromosomeElement* >::const_iterator 
            const_ce_iter;
        /// Chromosome reverse iterator
        typedef std::list< ChromosomeElement* >::reverse_iterator ce_riter;
        /// Chromosome const reverse iterator
        typedef std::list< ChromosomeElement* >::const_reverse_iterator 
            const_ce_riter;

        public:
        /// Constructor
        Chromosome();
        /// Constructor, the list is not copied. Memory responsabilities
        /// are moved to \c this.
        Chromosome( Genome *, std::list< ChromosomeElement* > * );
        /// Copy constructor
        explicit Chromosome( const Chromosome & );
        /// Destructor
        virtual ~Chromosome();
        /// Copy a chromosome into \c this
        void copy( const Chromosome & );
        /// Clone a chromosome
        Chromosome* clone() const;
        /// Back to the pool
        virtual bool toPool();

        /// Initialise (add interactions after reading in file)
        void initialise();
        /// All the mutations that can be handled within the chromosome
        /// are performed by invoking this method.
        int mutate();
        /// Copy and deletion events on the scale of genes.
        ce_iter geneMutate( ce_iter );
        
        /// Retrotransposons are copied around the genome (together with
        /// their accompanying LTRs) via reverse transcriptase. Removal of
        /// retroposons is done via reciprocal recombination.
        ce_iter retroposonMutate( ce_iter );
        /// Long Terminal Repeats are subject to double strand breaks.
        ce_iter repeatMutate( ce_iter );
        /// Interactions can indel or change weight or tag.
        ce_iter interactionMutate( ce_iter );
        /// Retrotransposons may arise newly at rare occasions.
        void retroposonInvade();
        /// Interactions are small binding sites that arise spontaneously
        void interactionGenesis();
        /// During lifetime genes may accidentally switch on or off
        int genePerturb();
        /// Select the upstream region of interactions of a certain gene
        ce_iter interactionSelect( ce_iter );
        
        /// Return an iterator to a random location on this chromosome.
        /// Retroposons are considered w/ LTRs, upstream regions may be split
        /// and, most importantly, end() is considered a location too.
        boost::tuple< Chromosome*, ce_iter > randChromosomeLocation();
        /// Return an iterator to a random element in the genome
        /// See \c randChromosomeElement
        boost::tuple< Chromosome*, ce_iter > randGenomeLocation();
        /// Return an iterator to a random gene in this chromosome
        boost::tuple< Chromosome*, ce_iter > randChromosomeGene();
        /// Return an iterator to a random gene in the genome
        boost::tuple< Chromosome*, ce_iter > randGenomeGene();
        /// Return all the segments
        std::list< Chromosome* > segments();
        /// Append a chromosome to the end of \c this
        void append( Chromosome* );
        
        /// Reset the state of all chromosome elements. It should be 
        /// performed before a new simulation update on a genome
        void reset();
        /// Empty the chromosome, the chromosome elements are not freed
        void clear();
        /// Make sure all flags of caching behaviour are set to update cache
        void recache();
        
        /// Set the parent (genome) of \c this
        void parent( Genome * );
        /// Get the parent of \c this
        Genome* parent() const;

        /// Get a reference to the contents of the chromosome
        const std::list< ChromosomeElement* > & elements() const;

        /// Set the state of the genes
        void setState( const boost::dynamic_bitset<> & );
        /// Set state of genes w/ tag
        void setState( int, int );
        
        /// Get nr retroposons
        uint nrRetroposons() const;
        /// Get nr of repeats
        uint nrRepeats() const;
        /// Get nr of dsbs
        uint nrDoubleStrandBreaks() const;
        /// Get nr of gene cp/rm
        std::vector< uint > nrMutations() const;
        /// Get the size of the chromosome (number of elements)
        uint size() const;
        /// Is the chromosome empty?
        bool empty() const;
        
        /// Write a (xml) representation to string
        virtual void writeXml( std::ostream & ) const;

        public:
        /// Zero all rates
        static void zeroRates();

        /// Set copy rate of a retroposon
        static void copyRetroposonRate( double );
        /// Get copy rate of a retroposon
        static double copyRetroposonRate();
        /// Set removal rate of a retroposon
        static void removeRetroposonRate( double );
        /// Get removal rate of a retroposon
        static double removeRetroposonRate();
        /// Set rate of repeat removal
        static void removeRepeatRate( double );
        /// Get rate of repeat removal
        static double removeRepeatRate();
        /// Set DSB recombination rate
        static void recombinationRate( double );
        /// Get DSB recombination rate
        static double recombinationRate();
        /// Set new retroposon rate
        static void newRetroposonRate( double );
        /// Get the new retroposon rate
        static double newRetroposonRate();
        
        /// Set gene copy rate
        static void copyGeneRate( double );
        /// Get gene copy rate
        static double copyGeneRate();
        /// Set gene deletion rate
        static void removeGeneRate( double );
        /// Get gene deletion rate
        static double removeGeneRate();
        /// Set tag rate
        static void tagGeneRate( double );
        /// Set threshold rate
        static void thresholdGeneRate( double );
        /// Set state switching rate (noise)
        static void stateRate( double );

        /// Set interaction weight toggle rate
        static void weightInteractionRate( double );
        /// Set interaction copy rate
        static void copyInteractionRate( double );
        /// Set interaction remove rate
        static void removeInteractionRate( double );
        /// Set interaction spontaneous genesis
        static void newInteractionRate( double );
        /// Set interaction tag change rate (to a random gene's tag)
        static void tagInteractionRate( double );
        
        public:
        class IsRetroposon : 
            public std::unary_function< ChromosomeElement*, bool > {
            public:
            IsRetroposon() {}
            bool operator()( ChromosomeElement *ce ) const
            { return typeid( *ce ) == typeid( Retroposon ); }
        };
        
        class IsGene : 
            public std::unary_function< ChromosomeElement*, bool > {
            public:
            IsGene() {}
            bool operator()( ChromosomeElement *ce ) const
            { return typeid( *ce ) == typeid( Gene ); }
        };
        
        /*class IsTransFactor :
            public std::unary_function< ChromosomeElement*, bool > {
            public:
            IsTransFactor() {}
            bool operator()( ChromosomeElement *ce ) const
            { return typeid( *ce ) == typeid( TransFactor ); }
        };
        
        class IsTargetGene :
            public std::unary_function< ChromosomeElement*, bool > {
            public:
            IsTargetGene() {}
            bool operator()( ChromosomeElement *ce ) const
            { return typeid( *ce ) == typeid( Gene ) and 
                typeid( *ce ) != typeid( TransFactor ); }
        };*/
        
        class IsInteraction : 
            public std::unary_function< ChromosomeElement*, bool > {
            public:
            IsInteraction() {}
            bool operator()( ChromosomeElement *ce ) const
            { return typeid( *ce ) == typeid( Interaction ); }
        };
        
        class IsRepeat :
            public std::unary_function< ChromosomeElement*, bool > {
            public:
            IsRepeat() {}
            bool operator()( ChromosomeElement *ce ) const
            { return typeid( *ce ) == typeid( Repeat ); }
        };
        
        class IsDoubleStrandBreak :
            public std::unary_function< ChromosomeElement*, bool > {
            public:
            IsDoubleStrandBreak() {}
            bool operator()( ChromosomeElement *ce ) const {
                Repeat *aux = dynamic_cast< Repeat* >( ce );
                if( aux ) {
                    return aux->hasDSB();
                } else {
                    return false;
                }
            }
        };

        private:
        // overloading list methods coz of length caching; insert, splice
        ce_iter insert( ce_iter, ChromosomeElement* );
        void splice( ce_iter, std::list< ChromosomeElement* >, uint, uint );
        
        private:
        Genome *parent_;
        std::list< ChromosomeElement* > *chro_;
        std::vector< uint > mut_events_;
        mutable uint retros_, ltrs_, len_;
        mutable bool update_retro_, update_ltr_, update_len_;
            
        // mutation rates 
        static double cp_tp_rate_;
        static double rm_tp_rate_;
        static double rm_ltr_rate_;
        static double new_tp_rate_;

        static double cp_gene_rate_;
        static double rm_gene_rate_;
        static double thr_gene_rate_;
        static double tag_gene_rate_;
        static double state_rate_;
        
        // rate of a dsb occurrence * rate of recombination repair
        static double dsb_recombination_;

        // gene mutation rates
        static double weight_ia_rate_;
        static double cp_ia_rate_;
        static double rm_ia_rate_;
        static double new_ia_rate_;
        static double tag_ia_rate_;
        
        public:
        // how many different mutations have we got?
        static int NR_MUT;
    };

    inline void Chromosome::parent( Genome *p )
    { parent_ = p; }

    inline Genome* Chromosome::parent() const 
    { return parent_; }

    inline boost::tuple< dorsalfin::Chromosome*, dorsalfin::Chromosome::ce_iter >
    Chromosome::randGenomeLocation()
    { return parent_->randLocation(); }

    inline boost::tuple< dorsalfin::Chromosome*, dorsalfin::Chromosome::ce_iter >
    Chromosome::randGenomeGene()
    { return parent_->randGene(); }
    
    inline const std::list< ChromosomeElement* >& Chromosome::elements() const
    { return *chro_; }

    inline uint Chromosome::nrDoubleStrandBreaks() const 
    { return mut_events_[ DSB ]; }

    inline std::vector< uint > Chromosome::nrMutations() const
    { return mut_events_; }
    
    inline bool Chromosome::empty() const 
    { return chro_->empty(); }
    
    inline void Chromosome::clear()
    { chro_->clear(); }
}
#endif

