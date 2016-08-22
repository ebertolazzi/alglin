/*--------------------------------------------------------------------------*\
 |                                                                          |
 |  Copyright (C) 2003                                                      |
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

#ifndef LU_BABD_CR_HH
#define LU_BABD_CR_HH

#include "Alglin.hh"
#include <vector>

// Eigen3
#include <Eigen/Dense>

#ifdef ALGLIN_USE_CXX11
  //#define LU_BABD_CR_USE_THREAD
#endif

#ifdef LU_BABD_CR_USE_THREAD
  #include <thread>
  #include <mutex>
  #include <condition_variable>
  #include <atomic>
  #define LU_BABD_CR_MAX_THREAD 256
#endif

namespace alglin {

  using namespace std ;

  /*
  //      _                        _ _       _    _   _ 
  //     / \   _ __ ___   ___   __| (_) ___ | |  | | | |
  //    / _ \ | '_ ` _ \ / _ \ / _` | |/ _ \| |  | | | |
  //   / ___ \| | | | | | (_) | (_| | | (_) | |__| |_| |
  //  /_/   \_\_| |_| |_|\___/ \__,_|_|\___/|_____\___/ 
  */

  //! LU decomposition of a BABD matrix
  /*!
   * 
   * \date     May 30, 2006
   * \version  1.0
   * \note     first release May 30, 2006
   *
   * \author   Enrico Bertolazzi
   *
   * \par      Affiliation:
   *           Department of Industrial Engineering<br>
   *           University of Trento<br>
   *           Via Sommarive 9, I-38123 Povo, Trento, Italy<br>
   *           `enrico.bertolazzi@unitn.it`
   *
   */
  template <typename t_Value, integer N>
  class CyclicReductionLU {
  private:

    typedef t_Value         valueType ;
    typedef t_Value*        valuePointer ;
    typedef t_Value const * valueConstPointer ;

    typedef Eigen::Matrix<integer,N,1> ivec_t ;
    typedef Eigen::Matrix<t_Value,N,1> vec_t ;
    typedef Eigen::Matrix<t_Value,N,N> mat_t ;
    typedef Eigen::Matrix<t_Value,Eigen::Dynamic,1>              dvec_t ;
    typedef Eigen::Matrix<t_Value,Eigen::Dynamic,Eigen::Dynamic> dmat_t ;
    //typedef Eigen::Block<mat_t>                                  mblock_t ;
    //typedef Eigen::Block<vec_t>                                  vblock_t ;
    //typedef Eigen::Block<ivec_t>                                 ivblock_t ;

    CyclicReductionLU(CyclicReductionLU const &) ;
    CyclicReductionLU const & operator = (CyclicReductionLU const &) ;

    integer nblock ; //!< total number of blocks
    integer q      ; //!< number final extra rows (q>=0)

    /*
    //
    //  Matrix structure
    //
    //                  n * nblock
    //    ________________^_________________
    //   /                                  \
    //     n     n     n                        n     q
    //  +-----+-----+-----+----.........-----+-----+-----+    \
    //  |  Ad | Au  |  0  |                  |  0  |  0  | n   |
    //  +-----+-----+-----+             -----+-----+-----+     |
    //  |  0  | Ad  | Au  |  0               |  0  |  0  | n   |
    //  +-----+-----+-----+-----+       -----+-----+-----+     |
    //  |  0  |  0  | Ad  | Au  |            |  0  |  0  | n   |
    //  +-----+-----+-----+-----+       -----+-----+-----+     |   
    //  |                                                :     |
    //  :                                                :      > n * nblock
    //  :                                                :     | 
    //  :                                                :     |
    //  :                                                :     |
    //  :                              +-----+-----+-----+     |
    //  :                              | Au  |  0  |  0  |     |
    //  :                        +-----+-----+-----+-----+     |
    //  :                        |  0  | Ad  | Au  |  0  | n   |
    //  +-----+-----+---......---+-----+-----+=====+=====+    /
    //  |     |     |            |     |     !     |     !
    //  | H0  |  0  |            |     |  0  ! HN  | Hq  ! m
    //  |     |     |            |     |     !     |     !
    //  +-----+-----+---......---+-----+-----+=====+=====+
    //
    */

    ///////////////////////////////////////////////////////

    std::vector<std::vector<bool> > LU_rows_blk ;
    std::vector<mat_t>              AdAuG_blk ;
    std::vector<ivec_t>             ipiv_blk ;
    mat_t EE, FF ;
    Eigen::FullPivLU<dmat_t>        LU_last ;

    mutable dvec_t tmpV ;
    mutable dvec_t tmpV1 ;

    integer NB ; // blocking factor

    #ifdef LU_BABD_CR_USE_THREAD
    mutable mutex              mtx0, mtx1, mtx2 ;
    mutable condition_variable cond0 ;
    mutable std::thread        threads[LU_BABD_CR_MAX_THREAD] ;
    mutable integer            to_be_done ;
    integer const              numThread ;
    #endif

    mutable integer k_block ;
    mutable integer jump_block ;

    integer
    LU_2_block( mat_t & A, mat_t & B, ivec_t & ipiv ) const ;

    #ifdef LU_BABD_CR_USE_THREAD
    void forward_reduce_mt( integer num_thread, vblock_t & y ) const ;
    void back_substitute_mt( integer num_thread, vblock_t & y ) const ;
    void reduction_mt( integer num_thread, integer nth ) ;
    #endif

    void forward_reduce( valuePointer y ) const ;
    void back_substitute( valuePointer y ) const ;
    void reduction() ;

  public:

    #ifdef LU_BABD_CR_USE_THREAD
    explicit CyclicReductionLU( integer nth = std::thread::hardware_concurrency() ) ;
    #else
    explicit CyclicReductionLU() ;
    #endif

    ~CyclicReductionLU() ;

    //! load matrix in the class
    /*!
      \param nblk number of (square) blocks
      \param n    size of the blocks
      \param q    extra bc
      \param AdAu pointer to the blocks diagonal ad upper diagonal
      \param H0   pointer to the block \f$ H_0 \f$
      \param HN   pointer to the block \f$ H_N \f$
      \param Hq   pointer to the block \f$ H_q \f$
      \code
      Matrix structure
      
                      n * nblock
        ________________^_________________
       /                                  \
         n     n     n                        n     q
      +-----+-----+-----+----.........-----+-----+-----+    \
      |  Ad | Au  |  0  |                  |  0  |  0  | n   |
      +-----+-----+-----+             -----+-----+-----+     |
      |  0  | Ad  | Au  |  0               |  0  |  0  | n   |
      +-----+-----+-----+-----+       -----+-----+-----+     |
      |  0  |  0  | Ad  | Au  |            |  0  |  0  | n   |
      +-----+-----+-----+-----+       -----+-----+-----+     |   
      |                                                :     |
      :                                                :      > n * nblock
      :                                                :     | 
      :                                                :     |
      :                                                :     |
      :                              +-----+-----+-----+     |
      :                              | Au  |  0  |  0  |     |
      :                        +-----+-----+-----+-----+     |
      :                        |  0  | Ad  | Au  |  0  | n   |
      +-----+-----+---......---+-----+-----+=====+=====+    /
      |     |     |            |     |     !     |     !
      | H0  |  0  |            |     |  0  ! HN  | Hq  ! m
      |     |     |            |     |     !     |     !
      +-----+-----+---......---+-----+-----+=====+=====+
      \endcode
    */
    void
    factorize( integer           nblk,
               integer           q,
               valueConstPointer AdAu,
               valueConstPointer H0,
               valueConstPointer HN,
               valueConstPointer Hq ) ;

    //! solve linear sistem using internal factorized matrix
    void
    solve( valuePointer in_out ) const ;

  } ;
}

#endif
