//
// Genome of an agent.
//
// by Anton Crombach, A.B.M.Crombach at uu.nl
//

#ifndef _DORSALFIN_GENOME_H_
#define _DORSALFIN_GENOME_H_

#include "defs.hh"

namespace dorsalfin {

    /// \class Genome
    /// \brief Container of chromosomes.
    ///
    /// The genome is one of the three levels in a cell that we distinguish. It 
    /// harbours a high level object model of a genome. The genome consists of
    /// chromosome (which consist of genes etc) and has several mutational 
    /// operators defined on them. On chromosome level the mutational process
    /// is chromosomal tail swapping. Such mutations may occur if double
    /// stranded breaks are not repaired correctly.
    class Genome : public XmlWriter {
        public:
        /// Chromosome iterator 
        typedef std::list< Chromosome* >::iterator chr_iter;
        /// Constant chromosome iterator
        typedef std::list< Chromosome* >::const_iterator const_chr_iter;

        public:
        /// (Dummy) constructor
        Genome();
        /// Constructor with list of chromosomes
        Genome( std::list< Chromosome* > * );
        /// Copy constructor
        explicit Genome( const Genome & );
        /// Get clone of this genome (OO pattern \a prototype)
        Genome* clone() const;
        /// Copy given genome into \c this
        void copy( const Genome & );
        /// Destructor
        virtual ~Genome();

        /// Initialise interactions according to network
        void initialise();
        /// Duplicate the genome
        void duplicate();
        /// Split the genome in two
        Genome* split();
        /// Mutate the genome
        int mutate();
        /// Perturb the network state
        void perturb();
        /// Write the genome to an output stream (in xml format)
        virtual void writeXml( std::ostream & ) const;

        /// Return a pointer to any location (respecting retroposons with
        /// accompanying LTRs)
        boost::tuple< dorsalfin::Chromosome*, 
            std::list< ChromosomeElement* >::iterator > randLocation();
        /// Return a pointer to a gene
        boost::tuple< dorsalfin::Chromosome*, 
            std::list< ChromosomeElement* >::iterator > randGene();
        
        /// Set state of the genes according to integer encoded network state
        void setState( const boost::dynamic_bitset<> & );
        /// Set state of genes with given tag
        void setState( int, int );
        /// Set states of genes with given tag
        void setState( const std::set< int > &, 
            const boost::dynamic_bitset<> & );
            
        /// Get number of chromosomes
        uint size() const;
        /// Get total number of chromosome elements (i.e. binding sites,
        /// transposons, module genes, ordinary genes)
        uint length() const;
        /// Is the genome empty?
        bool empty() const;

        /// Get constant reference to the chromosomes. Used by hopgraph.
        const std::list< Chromosome* > & chromosomes() const;
        /// Get nr of transposons
        uint nrRetroposons() const;
        /// Get nr of single ltrs
        uint nrRepeats() const;
        /// Get nr of double strand breaks
        uint nrDoubleStrandBreaks() const;
        /// Get nr of dsbs, gene cp/rm
        std::vector< uint > nrMutations() const;
        /// Any mutations at all?
        bool hasMutation() const;
        /// Get nr of perturbations
        uint nrPerturbations() const;

        private:
        void clear();
            
        private:
        std::list< Chromosome* > *chromos_;
        // The idea of the vector is that it is not changed during
        // the life of an individual, and only gets new values at mutating
        // (=reproduction)
        std::vector< uint > nr_mutations_;
        // Nr of perturbations
        uint nr_perturbed_;
    };
    
    inline uint Genome::size() const
    { return chromos_->size(); }

    inline bool Genome::empty() const
    { return chromos_->empty(); }

    inline const std::list< dorsalfin::Chromosome* > & 
    Genome::chromosomes() const 
    { return *chromos_; }
    
    inline uint Genome::nrPerturbations() const
    { return nr_perturbed_; }
}
#endif

