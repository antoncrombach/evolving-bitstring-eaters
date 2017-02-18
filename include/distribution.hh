//
// A few statistical distributions
//
// by Anton Crombach, A.B.M.Crombach at uu.nl
//

#ifndef _DORSALFIN_DISTRIBUTION_
#define _DORSALFIN_DISTRIBUTION_

namespace dorsalfin{

    double gammln( double xx );
    
    /// \class Weibull 
    /// \brief Weibull distribution. 
    ///
    /// A continuous probability distribution with a probability density 
    /// function 
    /// \f$\alpha \beta^{-\alpha} x^{\alpha - 1} e^{(x/\beta)^{\alpha}}\f$.
    ///
    /// Weibull distributions are usually used for the lifetime of objects. 
    /// In \c fluke the distribution is used for modelling the probability of 
    /// two transposons forming a composite transposon. This probability 
    /// decreases as the distance on the chromosome between the two transposons 
    /// increases. (obsolete)
    template< class T >
    class Weibull : public std::unary_function< T, T > {
        public:
            /// \arg alpha, shape parameter, float or double
            /// \arg beta, scale parameter, float or double
            Weibull( T alpha, T beta ) : alpha_( alpha ), beta_( beta ) {}
            // compiler generated cpy constr is ok with us

            /// Gets \c alpha.
            T alpha() { return alpha_; }
            /// Gets \c beta.
            T beta() { return beta_; }

            /// Returns probability at \c x.
            T operator()( T x ) const {
                return alpha_ * std::pow( beta_, -alpha_ ) * 
                    std::pow( x, alpha_ - 1 ) * 
                    std::exp( std::pow( -( x / beta_ ), alpha_ ) );
            }

        private:
            T alpha_, beta_;        
    };
    
    /// \class LinearDecrease
    /// \brief Linear decreasing probabilty distribution
    template< class T >
    class LinearDecrease : public std::unary_function< T, T > {
        public:
            /// \arg alpha, y value at x = 0
            /// \arg beta, x value at y = 0
            LinearDecrease( T alpha, T beta ) : alpha_( alpha ), beta_( beta ){}
            
            /// Gets \c alpha
            T alpha() { return alpha_; }
            /// Gets \c beta.
            T beta() { return beta_; }

            /// Returns probability at \c x.
            T operator()( T x ) const {
                return alpha_ - x * ( alpha_ / beta_ );
            }
            
        private:
            T alpha_, beta_;
    };

    template< typename RandomNumberGenerator, 
        typename T=int, typename F=double >
    class Binomial {
        public:
        Binomial() : n_old_( -1 ), p_old_( -1.0 ) {}
        // compiler generated cpy constr is ok with us
        
        /// Return a random deviate from a binomial distr with n trials and 
        /// probability of success pp.
        T operator()( T n, F pp, RandomNumberGenerator &unif ) {
            F am, em, p, bnl, sq, t, y;
            
            p = ( pp <= 0.5 ? pp: 1.0 - pp );
            am = n * p;
            if( n < 25 ) {
                // direct computation
                bnl = 0.0;
                for( int j = 1; j <= n; ++j ) {
                    if( unif() < p ) ++bnl;
                }
            } else if( am < 1.0 ) {
                // exponential / poisson
                int j;
                F g = exp( -am );
                F t = 1.0;
                for( j = 0; j <= n; ++j ) {
                    t *= unif();
                    if( t < g ) break;
                }
                bnl = ( j <= n ? j: n );
            } else {
                if( n != n_old_ ) {
                    en_ = n;
                    oldg_ = gammln( en_ + 1.0 );
                    n_old_ = n;
                }
                if( p != p_old_ ) {
                    pc_ = 1.0 - p;
                    plog_ = log( p );
                    pclog_ = log( pc_ );
                    p_old_ = p;
                }
                sq = sqrt( 2.0 * am * pc_ );
                do {
                    do {
                        y = tan( PI * unif() );
                        em = sq * y * am;
                    } while( em < 0.0 or em >= ( en_ + 1.0 ) );
                    em = floor( em );
                    t = 1.2 * sq * ( 1.0 + y * y ) * 
                        exp( oldg_ - gammln( em + 1.0 )
                        - gammln( en_ - em + 1.0 ) + em * plog_ +
                        ( en_ - em ) * pclog_ );
                } while( unif() > t );
                bnl = em;
            }
            if( p != pp ) bnl = n - bnl;
            return static_cast< int >( bnl + 0.5 );
        }
            
        private:
        T n_old_;
        F p_old_, pc_, plog_, pclog_, en_, oldg_;
    };

    template< typename T=double >
    class Gaussian : public std::unary_function< T, T > {
        public:
        Gaussian() : iset_( false ), gset_( 0.0 ) {}

        T operator()( T avg, T sd ) {
            T fac, rsq, v1, v2;
            
            if( !iset_ ) {
                do {
                    v1= 2.0 * uniform() - 1.0;
                    v2= 2.0 * uniform() - 1.0;
                    rsq = v1*v1 + v2*v2;
                } while( rsq >= 1.0 or rsq == 0.0 );
                fac = sqrt( -2.0 * log( rsq ) / rsq );
                gset_ = v1 * fac;
                iset_ = true;
                return sd * ( v2 * fac ) + avg;
            } else {
                iset_ = false;
                return sd * gset_ + avg;
            }
        }
        
        private:
        bool iset_;
        T gset_;
    };
}
#endif
