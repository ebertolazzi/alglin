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

#ifndef ABD_AMODIO_N_HH
#define ABD_AMODIO_N_HH

#include "Alglin.hh"
#include "Alglin++.hh"
#include <vector>

#ifdef ALGLIN_USE_CXX11
  #define BABD_AMODIO_N_USE_THREAD
#endif

#ifdef BABD_AMODIO_N_USE_THREAD
  #include <thread>
  #include <mutex>
  #include <condition_variable>
  #include <atomic>
  #define BABD_AMODIO_MAX_THREAD 256
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
  class AmodioN {
  private:
  
    typedef t_Value         valueType ;
    typedef t_Value*        valuePointer ;
    typedef t_Value const * valueConstPointer ;

    Malloc<valueType> baseValue ;
    Malloc<integer>   baseInteger ;
    std::vector<std::vector<bool> > LU_rows_blk ;

    AmodioN(AmodioN const &) ;
    AmodioN const & operator = (AmodioN const &) ;

    integer nblock ; //!< total number of blocks
    integer q      ; //!< number final rows (m>=n)
    integer m      ;
    integer nm     ;

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

    valuePointer G_blk ;
    valuePointer AdAu_blk ;
    valuePointer LU_blk ; // last LU and working space
    valuePointer tmpV ;
    valuePointer tmpM ;

    integer * ipiv_blk ;
    integer * LU_ipiv_blk ;

    mutable integer jump_block ;

    integer
    LU_2_block( valuePointer A,
                valuePointer B,
                integer      ipiv[] ) const ;

    #ifdef BABD_AMODIO_N_USE_THREAD
    mutable mutex              mtx0, mtx1, mtx2 ;
    mutable condition_variable cond0 ;
    mutable std::thread        threads[BABD_AMODIO_MAX_THREAD] ;
    mutable integer            to_be_done ;
            integer const      numThread ;
    mutable integer            usedThread ;
    mutable integer            jump_block_max_mt ;
    mutable valuePointer       y_thread ;

    void forward_reduce_mt( integer nth ) const ;
    void back_substitute_mt( integer nth ) const ;
    void reduction_mt( integer nth ) ;
    #endif

    void forward_reduce( valuePointer y ) const ;
    void back_substitute( valuePointer y, integer jump_block_min ) const ;
    void reduction() ;

  public:

    #ifdef BABD_AMODIO_N_USE_THREAD
    explicit AmodioN( integer nth = std::thread::hardware_concurrency() ) ;
    #else
    explicit AmodioN() ;
    #endif

    ~AmodioN() ;

    //! load matrix in the class
    /*!
      \param nblk number of (square) blocks
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
