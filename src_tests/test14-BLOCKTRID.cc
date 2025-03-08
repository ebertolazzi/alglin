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

#include "Alglin.hh"

#ifdef __clang__
#pragma clang diagnostic ignored "-Wexit-time-destructors"
#pragma clang diagnostic ignored "-Wglobal-constructors"
#endif

static Utils::Console msg( &std::cout,4);

using namespace std;

static
void
test1() {

  alglin::BlockTridiagonalSymmetic<alglin::real_type> BT;
  alglin::integer rBlocks[]{ 0, 2, 5, 7};
  BT.setup( 3, rBlocks );

  alglin::real_type D0[]{ 2, 1, 1, 1};
  alglin::real_type D1[]{ 2, 0, 1,
                          0, 0, 0,
                          1, 0, 2};
  alglin::real_type D2[]{ 2, 0,
                          0, 1};

  alglin::real_type L0[]{ 1, 1, 0,
                          0, 1, 1 };
  alglin::real_type L1[]{ 1, 1, 1,
                          1, 1, 1 };

  alglin::real_type rhs[]{
    5, 4, 6, 4, 6, 5, 4,
    5, 4, 6, 4, 6, 5, 4
  };

  BT.setD( 0, D0, 2 );
  BT.setD( 1, D1, 3 );
  BT.setD( 2, D2, 2 );
  BT.setL( 0, L0, 3 );
  BT.setL( 1, L1, 2 );

  BT.factorize( "BT" );

  //BT.solve( rhs );
  BT.solve( 2, rhs, 7 );

  for ( alglin::integer k{0}; k < rBlocks[3]; ++k )
    fmt::print("x[{}] = {}\n", k, rhs[7+k]);
}

static
void
test2() {

  alglin::BlockTridiagonalSymmetic<alglin::real_type> BT;
  alglin::integer rBlocks[]{ 0, 2, 5, 7};
  BT.setup( 3, rBlocks );

  alglin::integer ii[]{
    1, 2, 1, 2,
    3, 5, 3, 5,
    6, 7,
    3, 4, 4, 5,
    6, 6, 6, 7, 7, 7
  };

  alglin::integer jj[]{
    1, 1, 2, 2,
    3, 3, 5, 5,
    6, 7,
    1, 1, 2, 2,
    3, 4, 5, 3, 4, 5
  };

  alglin::real_type vals[]{
    2, 1, 1, 1,
    2, 1, 1, 2,
    2, 1,
    1, 1, 1, 1,
    1, 1, 1, 1, 1, 1
  };

  alglin::real_type rhs[]{
    5, 4, 6, 4, 6, 5, 4,
    5, 4, 6, 4, 6, 5, 4
  };

  BT.zero();
  for ( int k{0}; k < 20; ++k )
    BT(ii[k]-1,jj[k]-1) += vals[k];

  BT.factorize( "BT" );

  //BT.solve( rhs );
  BT.solve( 2, rhs, 7 );

  for ( alglin::integer k{0}; k < rBlocks[3]; ++k )
    fmt::print("x[{}] = {}\n", k, rhs[7+k]);
}


static
void
test3() {

  alglin::KKT<alglin::real_type> BT;

  alglin::integer rBlocks[]{ 0, 2, 4, 7 };
  alglin::integer ii[]{ 1,2,2,8,9,9,5,5,6,6,7,7,8,8,8,9,9,9,3,4,5,6,7,3,4,8,9,8,9 };

  alglin::integer jj[]{ 1,1,2,8,8,9,3,4,3,4,3,4,5,6,7,5,6,7,3,4,5,6,7,1,2,1,2,3,4 };

  alglin::real_type vals[]{2,1,1,1,1,1,1,0,2,0,2,0,5,0,0,0,6,0,2,2,3,3,3,-1,-1,1,2,3,4};

  alglin::real_type rhs[]{ -2, -9, -5, -10, 3, -21, 0, 28, 7 };

  BT.load_triblock( 7, 2, 3, rBlocks, vals, ii, -1, jj, -1, 29, true );
  //BT.factorize( 7, 2, 4, 4, vals, ii, -1, jj, -1, 29, true );

  //BT.zero();
  //for ( int k{0}; k < 29; ++k )
  //  BT(ii[k]-1,jj[k]-1) += vals[k];

  BT.factorize();
  BT.solve( rhs );

  for ( alglin::integer k{0}; k < rBlocks[3]; ++k )
    fmt::print("x[{}] = {}\n", k, rhs[k]);
}

#include <random>

static unsigned seed1 = 2;
static std::mt19937 generator(seed1);

static
alglin::real_type
rand( alglin::real_type xmin, alglin::real_type xmax ) {
  alglin::real_type random = alglin::real_type(generator())/generator.max();
  return xmin + (xmax-xmin)*random;
}

static
void
test4() {

  Utils::TicToc tm;

  alglin::KKT<alglin::real_type> BT1, BT2;

  std::vector<alglin::integer>   ii, jj;
  std::vector<alglin::real_type> vals, rhs1, rhs2;

  alglin::integer n      { 100000 };
  alglin::integer bksize { 50 };
  alglin::integer nblk   { n/bksize-4 };
  alglin::integer m      { n - bksize*nblk };

  for ( alglin::integer i{0}; i < n; ++i ) {
    ii.push_back( i );
    jj.push_back( i );
    vals.push_back( bksize );
    rhs1.push_back( i );
    rhs2.push_back( i );
    for ( alglin::integer j{-bksize}; j < 0; ++j ) {
      if ( i+j < 0 || i+j >= n ) continue;
      ii.push_back( i );
      jj.push_back( i+j );
      vals.push_back( rand(-1,1) );
    }
  }

  tm.tic();
  BT1.load_banded(
    n-m, m, bksize, bksize,
    &vals.front(),
    &ii.front(), 0,
    &jj.front(), 0,
    alglin::integer(vals.size()), true
  );
  tm.toc();
  fmt::print("LOAD1      = {:.5} [ms]\n",tm.elapsed_ms());

  tm.tic();
  BT2.load_triblock(
    n-m, m, nblk, bksize,
    &vals.front(),
    &ii.front(), 0,
    &jj.front(), 0,
    alglin::integer(vals.size()), true
  );

  tm.toc();
  fmt::print("LOAD2      = {:.5} [ms]\n", tm.elapsed_ms());

  tm.tic();
  BT1.factorize();
  tm.toc();
  fmt::print("factorize1 = {:.5} [ms]\n", tm.elapsed_ms());

  tm.tic();
  BT2.factorize();
  tm.toc();
  fmt::print("factorize2 = {:.5} [ms]\n", tm.elapsed_ms());

  tm.tic();
  BT1.solve( &rhs1.front() );
  tm.toc();
  fmt::print("solve1 = {:.5} [ms]\n", tm.elapsed_ms());

  tm.tic();
  BT2.solve( &rhs2.front() );
  tm.toc();
  fmt::print("solve2 = {:.5} [ms]\n", tm.elapsed_ms());

  alglin::real_type accerr = 0;
  alglin::real_type maxerr = 0;
  for ( size_t i{0}; i < size_t(n); ++i ) {
    alglin::real_type err = std::abs( rhs1[i] - rhs2[i] );
    accerr += err;
    if ( maxerr < err ) maxerr = err;
  }

  fmt::print("‖err‖₁ = {:.5} [ms]\n", accerr);
  fmt::print("‖err‖₁ = {:.5} [ms]\n", maxerr);

}

int
main() {
  cout << "test1\n";
  test1();
  cout << "\n\ntest2\n";
  test2();
  cout << "\n\ntest3\n";
  test3();
  cout << "\n\ntest4\n";
  test4();
  msg.green( "All done!\n" );
  return 0;
}
