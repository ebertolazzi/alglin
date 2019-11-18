/*--------------------------------------------------------------------------*\
 |                                                                          |
 |  Copyright (C) 2017                                                      |
 |                                                                          |
 |         , __                 , __                                        |
 |        /|/  \               /|/  \                                       |
 |         | __/ _   ,_         | __/ _   ,_                                |
 |         |   \|/  /  |  |   | |   \|/  /  |  |   |                        |
 |         |(__/|__/   |_/ \_/|/|(__/|__/   |_/ \_/|/                       |
 |                           /|                   /|                        |
 |                           \|                   \|                        |
 |                                                                          |
 |      Enrico Bertolazzi                                                   |
 |      Dipartimento di Ingegneria Industriale                              |
 |      Universita` degli Studi di Trento                                   |
 |      email: enrico.bertolazzi@unitn.it                                   |
 |                                                                          |
\*--------------------------------------------------------------------------*/

#ifndef BABD_BORDERED_CR_EIGEN3_HH
#define BABD_BORDERED_CR_EIGEN3_HH

#include "Alglin_Config.hh"
#include "Alglin_SuperLU.hh"
#include "Alglin_Eigen.hh"
#include "Alglin_threads.hh"

#include <iostream>

#ifndef BORDERED_CYCLIC_REDUCTION_MAX_THREAD
  #define BORDERED_CYCLIC_REDUCTION_MAX_THREAD 256
#endif

#ifdef __GNUC__ 
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wpadded"
#endif
#ifdef __clang__
#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wpadded"
#pragma clang diagnostic ignored "-Wc++98-compat"
#endif

namespace alglin {

  /*\
   |   ___             _                _    ___ ___
   |  | _ ) ___ _ _ __| |___ _ _ ___ __| |  / __| _ \
   |  | _ \/ _ \ '_/ _` / -_) '_/ -_) _` | | (__|   /
   |  |___/\___/_| \__,_\___|_| \___\__,_|  \___|_|_\
  \*/

  //! Cyclic reduction of a block bidiagonal matrix
  /*!
   *
   * \date     October 25, 2016
   * \version  1.0
   * \note     October 25, 2016
   *
   * \author   Enrico Bertolazzi
   *
   * \par      Affiliation:
   *           Department of Industrial Engineering<br>
   *           University of Trento <br>
   *           Via Sommarive 9, I-38123 Povo, Trento, Italy<br>
   *           enrico.bertolazzi\@unitn.it
   *
  \*/

  /*\
   |
   |  Matrix structure
   |
   |                 n * (nblock+1)
   |    ___________________^____________________
   |   /                                        \
   |    n   n   n                              n  qx  nx
   |  +---+---+---+----.................-----+---+---+---+   -+
   |  | D | E |   |                          |   |   | B | n  |
   |  +---+---+---+                     -----+---+---+---+    |
   |  |   | D | E |                          |   |   | B | n  |
   |  +---+---+---+---+                 -----+---+---+---+    |
   |  |   |   | D | E |                      |   |   | B | n  |
   |  +---+---+---+---+                 -----+---+---+---+    |
   |  :                                                  :    |
   |  :                                                  :    |
   |  :                                                  :     > n * nblock
   |  :                                                  :    |
   |  :                                                  :    |
   |  :                              +---+---+---+---+---+    |
   |  :                              | D | E |   |   | B | n  |
   |  :                              +---+---+---+---+---+    |
   |  :                                  | D | E |   | B | n  |
   |  +---+---+---................---+---+---+---+---+---+   -+
   |  |   |   |                          |   |   |   |   |    |
   |  |H0 | 0 |                          | 0 |HN | Hq| Hp|    | n+qr
   |  |   |   |                          |   |   |   |   |    |
   |  +---+---+---................---+---+---+---+---+---+   -+
   |  | C | C |                      | C | C | C | Cq| F |    | nr
   |  +---+---+---................---+---+---+---+---+---+   -+
   |                                             nr*qx
  \*/

  template <typename t_Value>
  class BorderedCR_eigen3 : public LinearSystemSolver<t_Value> {
  public:
    typedef t_Value valueType;

    typedef enum {
      BORDERED_LU      = 0,
      BORDERED_QR      = 1,
      BORDERED_QRP     = 3
    } BORDERED_Choice;

    typedef enum {
      BORDERED_LAST_LU   = 0,
      BORDERED_LAST_LUPQ = 1,
      BORDERED_LAST_QR   = 2,
      BORDERED_LAST_QRP  = 3,
      BORDERED_LAST_SVD  = 4,
      BORDERED_LAST_LSS  = 5,
      BORDERED_LAST_LSY  = 6
    } BORDERED_LAST_Choice;

  private:
    BorderedCR_eigen3(BorderedCR_eigen3 const &);
    BorderedCR_eigen3 const & operator = (BorderedCR_eigen3 const &);

  protected:

    Malloc<valueType> baseValue;
    Malloc<integer>   baseInteger;

    Malloc<valueType> superluValue;
    Malloc<int>       superluInteger;

    integer nblock; //!< total number of blocks
    integer n;      //!< size of square blocks
    integer qr, qx; //!< extra BC
    integer nr, nx; //!< border size

    // some derived constants
    integer n_x_2;
    integer n_x_n;
    integer n_x_nx;
    integer nr_x_n;
    integer nr_x_nx;
    integer Nr, Nc;
    integer Tsize;

    BORDERED_LAST_Choice last_selected;
    BORDERED_Choice      selected;

    void
    buildT(
      integer         nth,
      valueType const TOP[],
      valueType const BOTTOM[],
      valueType       T[],
      integer         iperm[]
    ) const;

    void
    applyT(
      integer         nth,
      valueType const T[],
      integer   const iperm[],
      valueType       TOP[],
      integer         ldTOP,
      valueType       BOTTOM[],
      integer         ldBOTTOM,
      integer         ncol
    ) const;

    void
    applyT(
      integer         nth,
      valueType const T[],
      integer   const iperm[],
      valueType       TOP[],
      valueType       BOTTOM[]
    ) const;

    // convert permutation to exchanges
    void
    permutation_to_exchange( integer nn, integer P[], integer S[] ) const {
      for ( integer i = 0; i < nn; ++i ) {
        integer j = i;
        while ( j < nn ) { if ( P[j] == i+1 ) break; ++j; }
        //LW_ASSERT0( j < nn, "permutation_to_exchange error!\n" );
        std::swap( P[j], P[i] );
        S[i] = j;
      }
    }

    /*
    //    __         _           _
    //   / _|__ _ __| |_ ___ _ _(_)______
    //  |  _/ _` / _|  _/ _ \ '_| |_ / -_)
    //  |_| \__,_\__|\__\___/_| |_/__\___|
    */

    void
    factorize_block( integer nth );

    void
    factorize_reduced();

    /*
    //    __                            _
    //   / _|___ _ ___ __ ____ _ _ _ __| |
    //  |  _/ _ \ '_\ V  V / _` | '_/ _` |
    //  |_| \___/_|  \_/\_/\__,_|_| \__,_|
    */

    void
    forward( integer nth, valueType x[], valueType xb[] ) const;

    void
    forward_n(
      integer   nth,
      integer   nrhs,
      valueType rhs[],
      integer   ldRhs
    ) const;

    void
    forward_reduced( valueType x[], valueType xb[] ) const;

    void
    forward_n_reduced(
      integer   nrhs,
      valueType rhs[],
      integer   ldRhs
    ) const;

    /*
    //   _             _                       _
    //  | |__  __ _ __| |____ __ ____ _ _ _ __| |
    //  | '_ \/ _` / _| / /\ V  V / _` | '_/ _` |
    //  |_.__/\__,_\__|_\_\ \_/\_/\__,_|_| \__,_|
    */

    void
    backward( integer nth, valueType x[] ) const;

    void
    backward_reduced( valueType x[] ) const;

    void
    backward_n(
      integer   nth,
      integer   nrhs,
      valueType rhs[],
      integer   ldRhs
    ) const;

    void
    backward_n_reduced(
      integer   nrhs,
      valueType rhs[],
      integer   ldRhs
    ) const;

    bool
    load_and_factorize_last();

    /*
    //   _         _
    //  | |__ _ __| |_
    //  | / _` (_-<  _|
    //  |_\__,_/__/\__|
    */

    bool
    solve_last( valueType [] ) const;

    bool
    solve_last(
      integer   nrhs,
      valueType rhs[],
      integer   ldRhs
    ) const;

    valueType * H0Nqp;
    valueType * Bmat;
    valueType * Cmat;
    valueType * Cqmat;
    valueType * Dmat;
    valueType * Emat;
    valueType * Fmat[BORDERED_CYCLIC_REDUCTION_MAX_THREAD];
    valueType * WorkT[BORDERED_CYCLIC_REDUCTION_MAX_THREAD];
    valueType * WorkQR[BORDERED_CYCLIC_REDUCTION_MAX_THREAD];

    // working block
    valueType * Tmat;
    valueType * Ttau;
    valueType * Work;
    integer   * Perm;
    integer     Lwork, LworkT, LworkQR;

    // last block
    valueType * Hmat;

    LU<valueType>   last_lu;
    LUPQ<valueType> last_lupq;
    QR<valueType>   last_qr;
    QRP<valueType>  last_qrp;
    SVD<valueType>  last_svd;
    LSS<valueType>  last_lss;
    LSY<valueType>  last_lsy;

    integer *iBlock, *kBlock;

    // used also with a unique thread
    integer maxThread, usedThread;
    mutable integer   * perm_thread[BORDERED_CYCLIC_REDUCTION_MAX_THREAD];
    mutable valueType * xb_thread[BORDERED_CYCLIC_REDUCTION_MAX_THREAD];
    mutable std::thread threads[BORDERED_CYCLIC_REDUCTION_MAX_THREAD];
    mutable SpinLock spin;

  public:

    using LinearSystemSolver<t_Value>::factorize;

    explicit
    BorderedCR_eigen3( integer nth = integer(1+std::thread::hardware_concurrency()/2) )
    : baseValue("BorderedCR_values")
    , baseInteger("BorderedCR_integers")
    , superluValue("BorderedCR_superluValue")
    , superluInteger("BorderedCR_superluInteger")
    , nblock(0)
    , n(0)
    , qr(0)
    , qx(0)
    , nr(0)
    , nx(0)
    , n_x_2(0)
    , n_x_n(0)
    , n_x_nx(0)
    , nr_x_n(0)
    , nr_x_nx(0)
    , Nr(0)
    , Nc(0)
    , last_selected(BORDERED_LAST_LU)
    , selected(BORDERED_LU)
    , H0Nqp(nullptr)
    , Bmat(nullptr)
    , Cmat(nullptr)
    , Cqmat(nullptr)
    , Dmat(nullptr)
    , Emat(nullptr)
    , Tmat(nullptr)
    , Ttau(nullptr)
    , Work(nullptr)
    , Perm(nullptr)
    , Lwork(0)
    , Hmat(nullptr)
    , maxThread(1)
    {
      LW_ASSERT(
        nth > 0 && nth <= BORDERED_CYCLIC_REDUCTION_MAX_THREAD,
        "Bad number of thread specification [{}]\n"
        "must be a number > 0 and <= {}\n",
        nth, BORDERED_CYCLIC_REDUCTION_MAX_THREAD
      );
      maxThread  = integer(std::thread::hardware_concurrency());
      usedThread = nth;
      #ifdef LAPACK_WRAPPER_USE_OPENBLAS
      openblas_set_num_threads(1);
      goto_set_num_threads(1);
      #endif
    }

    virtual
    ~BorderedCR_eigen3() ALGLIN_OVERRIDE
    {}

    void
    set_N_thread( integer nth ) {
      if      ( nth < 1         ) usedThread = 1;
      else if ( nth > maxThread ) usedThread = maxThread;
      else                        usedThread = nth;
    }

    integer get_used_thread() const { return usedThread; }
    integer get_max_thread()  const { return maxThread; }

    //! load matrix in the class
    void
    allocate(
      integer _nblock,
      integer _n,
      integer _qr,
      integer _qx,
      integer _nr,
      integer _nx
    );

    void
    dup( BorderedCR_eigen3 const & );

    /*!
     | \name Select Linar Algebra solver
     | @{
    \*/

    void select_LU()  { selected = BORDERED_LU; }
    void select_QR()  { selected = BORDERED_QR; }
    void select_QRP() { selected = BORDERED_QRP; }

    void select_last_LU()   { last_selected = BORDERED_LAST_LU;   }
    void select_last_LUPQ() { last_selected = BORDERED_LAST_LUPQ; }
    void select_last_QR()   { last_selected = BORDERED_LAST_QR;   }
    void select_last_QRP()  { last_selected = BORDERED_LAST_QRP;  }
    void select_last_SVD()  { last_selected = BORDERED_LAST_SVD;  }
    void select_last_LSS()  { last_selected = BORDERED_LAST_LSS;  }
    void select_last_LSY()  { last_selected = BORDERED_LAST_LSY;  }

    static
    std::string
    choice_to_string( BORDERED_Choice c ) {
      std::string res = "none";
      switch ( c ) {
      case BORDERED_LU:  res = "CyclicReduction+LU";  break;
      case BORDERED_QR:  res = "CyclicReduction+QR";  break;
      case BORDERED_QRP: res = "CyclicReduction+QRP"; break;
      }
      return res;
    }

    static
    std::string
    choice_to_string( BORDERED_LAST_Choice c ) {
      std::string res = "LastBlock not selected";
      switch ( c ) {
      case BORDERED_LAST_LU:   res = "LastBlock LU";   break;
      case BORDERED_LAST_LUPQ: res = "LastBlock LUPQ"; break;
      case BORDERED_LAST_QR:   res = "LastBlock QR";   break;
      case BORDERED_LAST_QRP:  res = "LastBlock QRP";  break;
      case BORDERED_LAST_SVD:  res = "LastBlock SVD";  break;
      case BORDERED_LAST_LSS:  res = "LastBlock LSS";  break;
      case BORDERED_LAST_LSY:  res = "LastBlock LSY";  break;
      }
      return res;
    }

    std::string
    info_algo() const {
      std::string a = choice_to_string(selected);
      std::string b = choice_to_string(last_selected);
      return a+" and "+b;
    }

    void info( ostream_type & stream ) const;

    /*!
     | @}
    \*/

    //! \brief Number of rows of the linear system
    integer
    numRows() const
    { return n * (nblock+1) + qx + nx; }

    //! \brief Number of columns of the linear system
    integer
    numCols() const
    { return n * (nblock+1) + qr + nr; }

    /*!
     | \name Filling all or part of the linear system with zero
     | @{
    \*/

    void zeroB()  { alglin::zero( nblock*n_x_nx, Cmat,    1 ); }
    void zeroD()  { alglin::zero( nblock*n_x_n,  Dmat,    1 ); }
    void zeroE()  { alglin::zero( nblock*n_x_n,  Emat,    1 ); }
    void zeroF()  { alglin::zero( nr_x_nx,       Fmat[0], 1 ); }
    void zeroH()  { alglin::zero( (n+qr)*Nc,     H0Nqp,   1 ); }
    void zeroC()  { alglin::zero( (nblock+1)*nr_x_n, Cmat, 1 ); }
    void zeroCq() { alglin::zero( nr*qx, Cqmat, 1 ); }

    void
    fillZero()
    { zeroB(); zeroC(); zeroCq(); zeroD(); zeroE(); zeroF(); zeroH(); }

    /*!
     |  @}
     |
     |  \name Access to single block
     |
     |  Matrix structure
     |
     |  \verbatim
     |
     |                 n * (nblock+1)
     |    ___________________^____________________
     |   /                                        \
     |    n   n   n                              n  qx  nx
     |  +---+---+---+----.................-----+---+---+---+   -+
     |  | D | E |   |                          |   |   | B | n  |
     |  +---+---+---+                     -----+---+---+---+    |
     |  |   | D | E |                          |   |   | B | n  |
     |  +---+---+---+---+                 -----+---+---+---+    |
     |  |   |   | D | E |                      |   |   | B | n  |
     |  +---+---+---+---+                 -----+---+---+---+    |
     |  :                                                  :    |
     |  :                                                  :    |
     |  :                                                  :     > n * nblock
     |  :                                                  :    |
     |  :                                                  :    |
     |  :                              +---+---+---+---+---+    |
     |  :                              | D | E |   |   | B | n  |
     |  :                              +---+---+---+---+---+    |
     |  :                                  | D | E |   | B | n  |
     |  +---+---+---................---+---+---+---+---+---+   -+
     |  |   |   |                          |   |   |   |   |    |
     |  |H0 | 0 |                          | 0 |HN | Hq| Hp|    | n+qr
     |  |   |   |                          |   |   |   |   |    |
     |  +---+---+---................---+---+---+---+---+---+   -+
     |  | C | C |                      | C | C | C | Cq| F |    | nr
     |  +---+---+---................---+---+---+---+---+---+   -+
     |                                             nr*qx
     |  \endverbatim
     |
     |  @{
    \*/

    // Border Right blocks
    /*\
     |  ____
     | | __ )
     | |  _ \
     | | |_) |
     | |____/
    \*/
    void
    loadB( integer nbl, valueType const B[], integer ldB )
    { gecopy( n, nx, B, ldB, Bmat + nbl*n_x_nx, n ); }

    void
    loadB( integer nbl, MatrixWrapper<valueType> const & B );

    void
    addtoB( integer nbl, valueType const B[], integer ldB ) {
      valueType * BB = Bmat + nbl*n_x_nx;
      geadd( n, nx, 1.0, B, ldB, 1.0, BB, n, BB, n );
    }

    void
    addtoB( integer nbl, MatrixWrapper<valueType> const & B );

    integer
    patternB( integer nbl, integer I[], integer J[], integer offs ) const;

    integer
    valuesB( integer nbl, valueType V[] ) const;

    // Border Bottom blocks
    /*\
     |   ____
     |  / ___|
     | | |
     | | |___
     |  \____|
    \*/
    void
    loadC( integer nbl, valueType const C[], integer ldC )
    { gecopy( nr, n, C, ldC, Cmat + nbl*nr_x_n, nr ); }

    void
    loadC( integer nbl, MatrixWrapper<valueType> const & C );

    void
    addtoC( integer nbl, valueType const C[], integer ldC ) {
      LW_ASSERT(
        ldC >= nr, "addtoC( {}, C, ldC = {} ) bad ldC\n", nbl, ldC
      );
      valueType * CC = Cmat + nbl*nr_x_n;
      geadd( nr, n, 1.0, C, ldC, 1.0, CC, nr, CC, nr );
    }

    void
    addtoC( integer nbl, MatrixWrapper<valueType> const & C );

    integer
    patternC( integer nbl, integer I[], integer J[], integer offs ) const;

    integer
    valuesC( integer nbl, valueType V[] ) const;

    // add to block nbl and nbl+1
    void
    addtoC2( integer nbl, valueType const C[], integer ldC ) {
      LW_ASSERT(
        ldC >= nr, "addtoC2( {}, C, ldC = {} ) bad ldC\n", nbl, ldC
      );
      valueType * CC = Cmat + nbl*nr_x_n;
      geadd( nr, n_x_2, 1.0, C, ldC, 1.0, CC, nr, CC, nr );
    }

    void
    addtoC2( integer nbl, MatrixWrapper<valueType> const & C );

    // -------------------------------------------------------------------------
    /*\
     |  ____
     | |  _ \
     | | | | |
     | | |_| |
     | |____/
    \*/
    void
    loadD( integer nbl, valueType const D[], integer ldD )
    { gecopy( n, n, D, ldD, Dmat + nbl*n_x_n, n ); }

    void
    loadD( integer nbl, MatrixWrapper<valueType> const & D );

    integer
    patternD( integer nbl, integer I[], integer J[], integer offs ) const;

    integer
    valuesD( integer nbl, valueType V[] ) const;

    /*\
     |  _____
     | | ____|
     | |  _|
     | | |___
     | |_____|
    \*/
    void
    loadE( integer nbl, valueType const E[], integer ldE )
    { gecopy( n, n, E, ldE, Emat + nbl*n_x_n, n ); }

    void
    loadE( integer nbl, MatrixWrapper<valueType> const & E );

    integer
    patternE( integer nbl, integer I[], integer J[], integer offs ) const;

    integer
    valuesE( integer nbl, valueType V[] ) const;

    void
    loadDE( integer nbl, valueType const DE[], integer ldDE ) {
      gecopy( n, n, DE, ldDE, Dmat + nbl*n_x_n, n ); DE += n*ldDE;
      gecopy( n, n, DE, ldDE, Emat + nbl*n_x_n, n );
    }

    void
    loadDEB( integer nbl, valueType const DEB[], integer ldDEB ) {
      gecopy( n, n,  DEB, ldDEB, Dmat + nbl*n_x_n,  n  ); DEB += n*ldDEB;
      gecopy( n, n,  DEB, ldDEB, Emat + nbl*n_x_n,  n  ); DEB += n*ldDEB;
      gecopy( n, nx, DEB, ldDEB, Bmat + nbl*n_x_nx, nx );
    }

    // -------------------------------------------------------------------------
    /*\
     |  _____
     | |  ___|
     | | |_
     | |  _|
     | |_|
    \*/
    void
    loadF( valueType const F[], integer ldF )
    { gecopy( nr, nx, F, ldF, Fmat[0], nr ); }

    void
    loadF( MatrixWrapper<valueType> const & F );

    void
    addtoF( valueType const F[], integer ldF )
    { gecopy( nr, nx, F, ldF, Fmat[0], nr ); }

    void
    addtoF( MatrixWrapper<valueType> const & F );

    integer
    patternF( integer I[], integer J[], integer offs ) const;

    integer
    valuesF( valueType V[] ) const;

    // -------------------------------------------------------------------------
    /*\
     |   ____
     |  / ___|__ _
     | | |   / _` |
     | | |__| (_| |
     |  \____\__, |
     |          |_|
    \*/
    void
    loadCq( valueType const Cq[], integer ldC )
    { gecopy( nr, qx, Cq, ldC, Cqmat, nr ); }

    void
    loadCq( MatrixWrapper<valueType> const & Cq );

    void
    loadCqF( valueType const CqF[], integer ldCF ) {
      gecopy( nr, qx, CqF, ldCF, Cqmat, nr ); CqF += qx*ldCF;
      gecopy( nr, nx, CqF, ldCF, Fmat[0],  nr );
    }

    integer
    patternCq( integer I[], integer J[], integer offs ) const;

    integer
    valuesCq( valueType V[] ) const;

    // -------------------------------------------------------------------------
    /*\
     |  _   _
     | | | | |
     | | |_| |
     | |  _  |
     | |_| |_|
    \*/
    integer
    patternH( integer I[], integer J[], integer offs ) const;

    integer
    valuesH( valueType V[] ) const;

    void
    loadBottom(
      valueType const H0[], integer ld0,
      valueType const HN[], integer ldN,
      valueType const Hq[], integer ldQ,
      valueType const Hp[], integer ldP
    );

    void
    loadBottom(
      MatrixWrapper<valueType> const & H0,
      MatrixWrapper<valueType> const & HN,
      MatrixWrapper<valueType> const & Hq,
      MatrixWrapper<valueType> const & Hp
    );

    void
    loadBottom( valueType const _H0Nqp[], integer ldH ) {
      integer nq = n+qr;
      gecopy( nq, Nc, _H0Nqp, ldH, H0Nqp, nq );
    }

    void
    loadBottom( MatrixWrapper<valueType> const & H );

    /*\
     |  +---+---+---+---+
     |  | C | C | Cq| F |
     |  +---+---+---+---+
    \*/

    void
    loadBottom2(
      valueType const C0[], integer ld0,
      valueType const CN[], integer ldN,
      valueType const Cq[], integer ldCq,
      valueType const F[],  integer ldF
    );

    void
    loadBottom2(
      MatrixWrapper<valueType> const & C0,
      MatrixWrapper<valueType> const & CN,
      MatrixWrapper<valueType> const & Cq,
      MatrixWrapper<valueType> const & F
    );

    void
    loadBottom2( MatrixWrapper<valueType> const & H );

    /*!
     | @}
    \*/

    /*!
     |
     | \name Access to blocks by element
     |
     | @{
    \*/

    t_Value &
    B( integer nbl, integer i, integer j )
    { return Bmat[ nbl*n_x_nx + i + j*n ]; }

    t_Value const &
    B( integer nbl, integer i, integer j ) const
    { return Bmat[ nbl*n_x_nx + i + j*n ]; }

    t_Value &
    C( integer nbl, integer i, integer j )
    { return Cmat[ nbl*nr_x_n + i + j*nr ]; }

    t_Value const &
    C( integer nbl, integer i, integer j ) const
    { return Cmat[ nbl*nr_x_n + i + j*nr ]; }

    t_Value &
    D( integer nbl, integer i, integer j )
    { return Dmat[ nbl*n_x_n + i + j*n]; }

    t_Value const &
    D( integer nbl, integer i, integer j ) const
    { return Dmat[ nbl*n_x_n + i + j*n]; }

    t_Value &
    E( integer nbl, integer i, integer j )
    { return Emat[ nbl*n_x_n + i + j*n]; }

    t_Value const &
    E( integer nbl, integer i, integer j ) const
    { return Emat[ nbl*n_x_n + i + j*n]; }

    t_Value &
    F( integer i, integer j )
    { return Fmat[0][ i + j*nr ]; }

    t_Value const &
    F( integer i, integer j ) const
    { return Fmat[0][ i + j*nr ]; }

    t_Value &
    Cq( integer i, integer j )
    { return Cqmat[ i + j*nr ]; }

    t_Value const &
    Cq( integer i, integer j ) const
    { return Cqmat[ i + j*nr ]; }

    t_Value &
    H( integer i, integer j )
    { return H0Nqp[ i + j*(n+qr) ]; }

    t_Value const &
    H( integer i, integer j ) const
    { return H0Nqp[ i + j*(n+qr) ]; }

    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    void
    B( integer nbl, MatrixWrapper<valueType> & B_wrap )
    { B_wrap.setup( Bmat + nbl*n_x_nx, n, nx, n ); }

    void
    C( integer nbl, MatrixWrapper<valueType> & C_wrap )
    { C_wrap.setup( Cmat + nbl*nr_x_n, nr, n, nr); }

    void
    D( integer nbl, MatrixWrapper<valueType> & D_wrap )
    { D_wrap.setup( Dmat + nbl*n_x_n, n, n, n ); }

    void
    E( integer nbl, MatrixWrapper<valueType> & E_wrap )
    { E_wrap.setup( Emat + nbl*n_x_n, n, n, n ); }

    void
    F( MatrixWrapper<valueType> & F_wrap )
    { F_wrap.setup( Fmat[0], nr, nx, nr); }

    void
    Cq( MatrixWrapper<valueType> & Cq_wrap )
    { Cq_wrap.setup( Cqmat, nr, qx, nr ); }

    void
    H( MatrixWrapper<valueType> & H_wrap ) const
    { H_wrap.setup( H0Nqp, n + qr, Nc, n + qr ); }

    /*!
     | @}
    \*/

    void factorize();

    /*\
     |         _      _               _
     |  __   _(_)_ __| |_ _   _  __ _| |___
     |  \ \ / / | '__| __| | | |/ _` | / __|
     |   \ V /| | |  | |_| |_| | (_| | \__ \
     |    \_/ |_|_|   \__|\__,_|\__,_|_|___/
    \*/

    virtual
    bool
    solve( valueType x[] ) const ALGLIN_OVERRIDE;

    virtual
    bool
    solve( integer nrhs, valueType rhs[], integer ldRhs ) const ALGLIN_OVERRIDE;

    virtual
    bool
    t_solve( valueType [] ) const ALGLIN_OVERRIDE {
      LW_ERROR0( "BorderedCR::t_solve() not defined\n" );
    }

    virtual
    bool
    t_solve( integer, valueType [], integer ) const ALGLIN_OVERRIDE {
      LW_ERROR0( "BorderedCR::t_solve() not defined\n" );
    }

    /*\
     |     _
     |    / \  _   ___  __
     |   / _ \| | | \ \/ /
     |  / ___ \ |_| |>  <
     | /_/   \_\__,_/_/\_\
     |
    \*/
    void
    Mv( valueType const x[], valueType res[] ) const {
      alglin::zero( numRows(), res, 1 );
      addMv( x, res );
    }

    void
    addMv( valueType const x[], valueType res[] ) const;

    /*\
     |  ____
     | / ___| _ __   __ _ _ __ ___  ___
     | \___ \| '_ \ / _` | '__/ __|/ _ \
     |  ___) | |_) | (_| | |  \__ \  __/
     | |____/| .__/ \__,_|_|  |___/\___|
     |       |_|
     |
    \*/

    integer
    sparseNnz() const {
      return nblock*(2*n_x_n+n_x_nx+nr_x_n) +
             nr_x_n + nr*(qx+nx) + (n+qr)*(2*n+qx+nx);
    }

    void
    sparsePattern( integer I[], integer J[], integer offs ) const;

    void
    sparseValues( valueType V[] ) const;

    void
    sparseLoad(
      valueType const M_values[],
      integer   const M_row[], integer r_offs,
      integer   const M_col[], integer c_offs,
      integer         M_nnz
    );

  };

  // explicit instantiation declaration to suppress warnings

  #ifdef __clang__
  #pragma clang diagnostic push
  #pragma clang diagnostic ignored "-Wc++98-compat-pedantic"
  #pragma clang diagnostic ignored "-Wweak-template-vtables"
  #endif

  extern template class BorderedCR_eigen3<float>;
  extern template class BorderedCR_eigen3<double>;

  #ifdef __clang__
  #pragma clang diagnostic pop
  #endif
}

#endif
