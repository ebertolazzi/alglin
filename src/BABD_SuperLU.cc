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

#include "AlglinConfig.hh"

#ifdef __GNUC__ 
#pragma GCC diagnostic ignored "-Wsign-conversion"
#endif
#ifdef __clang__
#pragma clang diagnostic ignored "-Wsign-conversion"
#endif

#include "BABD_SuperLU.hh"

#include <iostream>
#include <algorithm>

using namespace std;

namespace alglin {

  using namespace std;

  BABD_SuperLU::BABD_SuperLU()
  : baseValue("SuperLU_values")
  , baseInteger("SuperLU_integes")
  { }

  BABD_SuperLU::~BABD_SuperLU() {
    baseValue   . free();
    baseInteger . free();
  }

  /*\
   |    __            _             _
   |   / _| __ _  ___| |_ ___  _ __(_)_______
   |  | |_ / _` |/ __| __/ _ \| '__| |_  / _ \
   |  |  _| (_| | (__| || (_) | |  | |/ /  __/
   |  |_|  \__,_|\___|\__\___/|_|  |_/___\___|
  \*/
  void
  BABD_SuperLU::factorize(
    integer         _nblock,
    integer         _n,
    integer         q,
    valueType const AdAu[],
    valueType const H0[],
    valueType const HN[],
    valueType const Hq[]
  ) {

    nblock = _nblock;
    n      = _n;
    m      = _n+q;
    nnz    = 2*nblock*n*n + (n+m)*m;
    neq    = nblock*n + m;

    baseInteger.allocate(size_t(2*nnz+4*neq+1));
    baseValue.allocate(size_t(nnz+2*neq));

    set_default_options(&slu_options);

    // Initialize the statistics variables.
    StatInit(&slu_stats);

    perm_r = baseInteger(size_t(neq)); /* row permutations from partial pivoting */
    perm_c = baseInteger(size_t(neq)); /* column permutation vector */
    etree  = baseInteger(size_t(neq));

    double * values = baseValue(size_t(nnz));
    double * sumR   = baseValue(size_t(neq));
    double * sumC   = baseValue(size_t(neq));
    int    * rowind = baseInteger(size_t(nnz));
    int    * colptr = baseInteger(size_t(neq+1));

    integer kk = 0;
    integer jj = 0;
    integer ee = nblock*n;
    colptr[0] = 0;
    for ( integer k = 0; k <= nblock; ++k ) {
      valueType const * Ad = AdAu + 2*n*n*k;
      valueType const * Au = Ad - n*n;
      integer ii = k*n;
      for ( integer j = 0; j < n; ++j ) {
        if ( k > 0       ) for ( integer i = 0; i < n; ++i ) { rowind[kk] = i+ii-n; values[kk] = Au[i+j*n]; ++kk; }
        if ( k < nblock  ) for ( integer i = 0; i < n; ++i ) { rowind[kk] = i+ii; values[kk] = Ad[i+j*n]; ++kk; }
        if ( k == 0      ) for ( integer i = 0; i < m; ++i ) { rowind[kk] = i+ee; values[kk] = H0[i+j*m]; ++kk; }
        if ( k == nblock ) for ( integer i = 0; i < m; ++i ) { rowind[kk] = i+ee; values[kk] = HN[i+j*m]; ++kk; }
        colptr[++jj] = kk;
      }
    }
    for ( integer j = 0; j < q; ++j ) {
      for ( integer i = 0; i < m; ++i ) { rowind[kk] = i+ee; values[kk++] = Hq[i+j*m]; }
      colptr[++jj] = kk;
    }

    ALGLIN_ASSERT( jj == neq, "SuperLU::factorize -- bad matrix format" );

    // Create matrix A in the format expected by SuperLU.
    //cout << "Create matrix A in the format expected by SuperLU.\n";
    dCreate_CompCol_Matrix( &A, neq, neq, nnz,
                            values, rowind, colptr,
                            SLU_NC, SLU_D, SLU_GE );

    std::fill( sumR, sumR+neq, 0 );
    std::fill( sumC, sumC+neq, 0 );
    for ( integer j = 0; j < neq; ++j ) {
      for ( integer ii = colptr[j]; ii < colptr[j+1]; ++ii ) {
        integer i = rowind[ii];
        double absAij = std::abs(values[ii]);
        sumR[i] += absAij;
        sumC[j] += absAij;
      }
    }
    one_norm_A = *std::max_element( sumC, sumC+neq );
    inf_norm_A = *std::max_element( sumR, sumR+neq );

    #if 0
    ofstream file("Jac.txt");
    file << "A := Matrix(" << neq << ", " << neq << ", fill = 0);\n";
    for ( integer j = 0; j < neq; ++j ) {
      for ( integer ii = colptr[j]; ii < colptr[j+1]; ++ii ) {
        integer i = rowind[ii];
        file << "A[" << i+1 << "," << j+1 << "] := " << values[ii] << ";\n";
      }
    }
    file.close();
    #endif

    factorize();
  }

  /*\
   |    __            _             _
   |   / _| __ _  ___| |_ ___  _ __(_)_______
   |  | |_ / _` |/ __| __/ _ \| '__| |_  / _ \
   |  |  _| (_| | (__| || (_) | |  | |/ /  __/
   |  |_|  \__,_|\___|\__\___/|_|  |_/___\___|
  \*/
  void
  BABD_SuperLU::factorize() {

    SuperMatrix AC;

    #if defined(SUPERLU_MAJOR_VERSION) && SUPERLU_MAJOR_VERSION >= 5
    GlobalLU_t  glu;
    #endif

    /*
     * Get column permutation vector perm_c[], according to permc_spec:
     *   ColPerm = 0: natural ordering
     *   ColPerm = 1: minimum degree on structure of A'*A
     *   ColPerm = 2: minimum degree on structure of A'+A
     *   ColPerm = 3: approximate minimum degree for unsymmetric matrices
     */
    //cout << "get_perm_c.\n";
    get_perm_c(slu_options.ColPerm, &A, perm_c);
    //cout << "sp_preorder.\n";
    sp_preorder(&slu_options, &A, perm_c, etree, &AC);

    int panel_size = sp_ienv(1);
    int relax      = sp_ienv(2);
    int info = 0;
    //cout << "dgstrf.\n";
    dgstrf(
      &slu_options, &AC, relax, panel_size,
      etree, nullptr, 0, perm_c, perm_r, &L, &U,
    #if defined(SUPERLU_MAJOR_VERSION) && SUPERLU_MAJOR_VERSION >= 5
      &glu,
    #endif
      &slu_stats, &info
    );

    // Free un-wanted storage
    Destroy_SuperMatrix_Store(&A);
    Destroy_CompCol_Permuted(&AC);
    StatFree(&slu_stats);

    ALGLIN_ASSERT( info == 0, "BABD_SuperLU::factorize -- dgstrf() error returns INFO= " << info );
    //cout << "done\n";
  }

  /*\
   |             _
   |   ___  ___ | |_   _____
   |  / __|/ _ \| \ \ / / _ \
   |  \__ \ (_) | |\ V /  __/
   |  |___/\___/|_| \_/ \___|
  \*/

  void
  BABD_SuperLU::solve( valueType y[] ) const {

    int         info;
    SuperMatrix B;
    trans_t     trans = NOTRANS;

    // Initialize the statistics variables.
    StatInit(&slu_stats);
    dCreate_Dense_Matrix(
      &B, L.nrow /* neq */,
      1 /* nrhs */,
      y,
      L.nrow /* ldy */,
      SLU_DN, SLU_D, SLU_GE
    );

    // Solve the system A*X=B, overwriting B with X.
    dgstrs( trans, &L, &U, perm_c, perm_r, &B, &slu_stats, &info);

    Destroy_SuperMatrix_Store(&B);
    StatFree(&slu_stats);

    ALGLIN_ASSERT( info == 0, "BABD_SuperLU::solve -- dgstrs() error returns INFO= " << info );
  }

  /*\
   |                       _
   |    ___ ___  _ __   __| |
   |   / __/ _ \| '_ \ / _` |
   |  | (_| (_) | | | | (_| |
   |   \___\___/|_| |_|\__,_|
  \*/
  void
  BABD_SuperLU::cond( valueType & rcond_1, valueType & rcond_inf ) const {
    int info1 = 0, info2 = 0;
    StatInit(&slu_stats);
    dgscon(	const_cast<char*>("1"), &L, &U, one_norm_A, &rcond_1, &slu_stats, &info1 );
    if ( info1 == 0 ) dgscon( const_cast<char*>("I"), &L, &U, inf_norm_A, &rcond_inf, &slu_stats, &info2 );
    StatFree(&slu_stats);
    ALGLIN_ASSERT( info1 == 0, "BABD_SuperLU::cond -- dgscon() error returns INFO= " << info1 );
    ALGLIN_ASSERT( info2 == 0, "BABD_SuperLU::cond -- dgscon() error returns INFO= " << info2 );
  }

}
