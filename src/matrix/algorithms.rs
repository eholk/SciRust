use std::num;

use matrix::{BasicMatrix, Create, SubMatrix, 
             SubMatrix_t, TransposeMatrix, Vector,
             col, row};
use matrix::generate::zero_matrix;
use matrix::util::tracefn;

//pub mod par;

pub fn dot<T: num::Num, L: Vector<T>, R: Vector<T>>(lhs: &L, rhs: &R) -> T {
    do tracefn(("dot", lhs.len(), rhs.len())) {
    assert!(lhs.len() > 0)
    if lhs.len() != rhs.len() {
        fail!(~"Invalid vector lengths.")
    }

    //error!("%? ### %?", lhs, rhs);

    info!("a");
    let a = lhs[0];
    info!("b");
    let b = rhs[0];

    let mut acc : T = a * b;
    for i in range(1, lhs.len()) {
        acc = acc + (lhs[i] * rhs[i])
    }

    acc
    }
}

pub fn mat_mul<T: Num + num::FromPrimitive, LHS: BasicMatrix<T>, RHS: BasicMatrix<T>, Res: BasicMatrix<T> + Create<T>> (lhs: &LHS, rhs: &RHS) -> Res
{
    do tracefn(("mat_mul",
               (lhs.num_rows(), lhs.num_cols()),
               (rhs.num_rows(), rhs.num_cols()))) {
    if lhs.num_cols() != rhs.num_rows() {
        fail!(format!("Incompatible matrix sizes. LHS: {:?}, RHS: {:?}",
                   (lhs.num_rows(), lhs.num_cols()),
                   (rhs.num_rows(), rhs.num_cols())))
    }

    Create::<T>::create(lhs.num_rows(), rhs.num_cols(),
                        |i, j| dot(&row(lhs, i), &col(rhs, j)))
    }
}

// M -> (A, B, C, D)
fn subdivide<'a, T: Num, N: BasicMatrix<T>,
             M: BasicMatrix<T> + SubMatrix<T, N>>
            (M: &'a M)
    -> (N, N, N, N)
{
    let H = M.num_rows();
    let W = M.num_cols();

    let H2  = H / 2;
    let H2a = H - H2;
    let W2  = W / 2;
    let W2a = W - W2;

    let A = SubMatrix(M,  0,  0,  H2,  W2);
    let B = SubMatrix(M,  0, W2,  H2, W2a);
    let C = SubMatrix(M, H2,  0, H2a,  W2);
    let D = SubMatrix(M, H2, W2, H2a, W2a);

    (A, B, C, D)
}

pub fn mat_mul_blocked
<T: Num + num::FromPrimitive,
SM: BasicMatrix<T>,
LHS: BasicMatrix<T> + SubMatrix<T, SM>, RHS: BasicMatrix<T> + SubMatrix<T, SM>,
Res: BasicMatrix<T> + Create<T> + SubMatrix<T, SM>>
(lhs: &LHS, rhs: &RHS)
 -> Res
{
    if lhs.num_cols() != rhs.num_rows() {
        fail!(format!("Incompatible matrix sizes. LHS: {:?}, RHS: {:?}",
                   (lhs.num_rows(), lhs.num_cols()),
                   (rhs.num_rows(), rhs.num_cols())))
    }

    static CUTOFF: uint = 32;
    
    if     lhs.num_rows() <= CUTOFF
        || lhs.num_cols() <= CUTOFF
        || rhs.num_rows() <= CUTOFF
        || rhs.num_cols() <= CUTOFF
    {
        mat_mul(lhs, rhs)
    }
    else {
        let (A, B, C, D) = subdivide(lhs);
        let (E, F, G, H) = subdivide(rhs);

        let A: Res = convert(&A);
        let B: Res = convert(&B);
        let C: Res = convert(&C);
        let D: Res = convert(&D);
        
        let E: Res = convert(&E);
        let F: Res = convert(&F);
        let G: Res = convert(&G);
        let H: Res = convert(&H);

        let AE: Res = mat_mul_blocked(&A, &E);
        let BG: Res = mat_mul_blocked(&B, &G);

        let An: Res = mat_add(&AE, &BG);

        let AF: Res = mat_mul_blocked(&A, &F);
        let BH: Res = mat_mul_blocked(&B, &H);

        let Bn: Res = mat_add(&AF, &BH);

        let CE: Res = mat_mul_blocked(&C, &E);
        let DG: Res = mat_mul_blocked(&D, &G);

        let Cn: Res = mat_add(&CE, &DG);

        let CF: Res = mat_mul_blocked(&C, &F);
        let DH: Res = mat_mul_blocked(&D, &H);

        let Dn: Res = mat_add(&CF, &DH);

        let top: Res = concat_cols(&An, &Bn);
        let bot: Res = concat_cols(&Cn, &Dn);

        concat_rows(&top, &bot)
    }
}

pub fn mat_add<T: Num, LHS: BasicMatrix<T>, RHS: BasicMatrix<T>, Res: BasicMatrix<T> + Create<T>> (lhs: &LHS, rhs: &RHS) -> Res
{
    if lhs.num_cols() != rhs.num_cols() || lhs.num_rows() != rhs.num_rows() {
        fail!(~"Incompatible matrix sizes")
    }

    do Create::<T>::create(lhs.num_rows(), rhs.num_cols()) |i, j| {
        lhs.get(i, j) + rhs.get(i, j)
    }
}

pub fn mat_sub<T: Num, LHS: BasicMatrix<T>, RHS: BasicMatrix<T>, Res: BasicMatrix<T> + Create<T>> (lhs: &LHS, rhs: &RHS) -> Res
{
    if lhs.num_cols() != rhs.num_cols() || lhs.num_rows() != rhs.num_rows() {
        fail!(format!("Incompatible matrix sizes. LHS: {:?}, RHS: {:?}",
                   (lhs.num_rows(), lhs.num_cols()),
                   (rhs.num_rows(), rhs.num_cols())))
    }

    do Create::<T>::create(lhs.num_rows(), rhs.num_cols()) |i, j| {
        lhs.get(i, j) - rhs.get(i, j)
    }
}

pub fn mat_x_inplace<T: Num, M: BasicMatrix<T>>(A: &mut M, x: T) {
    for_each(A, |_i, _j, y| x * y);
}

pub fn transpose<T, M: BasicMatrix<T>, R: BasicMatrix<T> + Create<T>>(m: &M) -> R {
    do Create::<T>::create(m.num_cols(), m.num_rows()) |i, j| {
        m.get(j, i)
    }
}

pub fn cholesky_seq_inplace_raw<M: BasicMatrix<f64>>(A: &mut M, start: uint) {
    assert!(A.num_rows() == A.num_cols());
    let N = A.num_rows();
    for k in range(start, N) {
        let Akk = A.get(k, k);

        A.set(k, k, Akk.sqrt());

        for i in range(k + 1, N) {
            let Aik = A.get(i, k);
            A.set(i, k, Aik / Akk);
        }

        for i in range(k + 1, N) {
            let Aik = A.get(i, k);
            for j in range(k + 1, i + 1) {
                let Ajk = A.get(j, k);
                let Aij = A.get(i, j);
                A.set(i, j, Aij - Aik * Ajk);
            }
        }
    }
}

pub fn for_each<T, M: BasicMatrix<T>>(A: &mut M, f: &fn(uint, uint, T) -> T) {
    for i in range(0, A.num_rows()) {
        for j in range(0, A.num_cols()) {
            let old = A.get(i, j);
            A.set(i, j, f(i, j, old))
        }
    }
}

pub fn cholesky_seq_inplace<M: BasicMatrix<f64>>(A: &mut M) {
    cholesky_seq_inplace_start::<M>(A, 0);
}

pub fn cholesky_seq_inplace_start<M: BasicMatrix<f64>>(A: &mut M,
                                                       start: uint) {
    cholesky_seq_inplace_raw::<M>(A, start);

    for i in range(start, A.num_rows()) {
        for j in range(i + 1, A.num_cols()) {
            A.set(i, j, num::from_int(0).unwrap())
        }
    }    
}


pub fn concat_rows<T, LHS: BasicMatrix<T>, RHS: BasicMatrix<T>, R: BasicMatrix<T> + Create<T>>(A: &LHS, B: &RHS) -> R {
    assert!(A.num_cols() == B.num_cols());

    let N = A.num_rows();
    do Create::<T>::create(N + B.num_rows(), A.num_cols()) |i, j| {
        if i < N {
            A.get(i, j)
        }
        else {
            B.get(i - N, j)
        }
    }
}

pub fn concat_cols<T, LHS: BasicMatrix<T>, RHS: BasicMatrix<T>, R: BasicMatrix<T> + Create<T>>(A: &LHS, B: &RHS) -> R {
    //error!("concat: {:?}, {:?}",
    //       (A.num_rows(), A.num_cols()),
    //       (B.num_rows(), B.num_cols()));
           
    assert!(A.num_rows() == B.num_rows());

    let N = A.num_cols();

    do Create::<T>::create(B.num_rows(), N + B.num_cols()) |i, j| {
        if j < N {
            A.get(i, j)
        }
        else {
            B.get(i, j - N)
        }
    }
}

pub fn convert<T, M: BasicMatrix<T>, R: BasicMatrix<T> + Create<T>>(M: &M) -> R {
    Create::<T>::create(M.num_rows(), M.num_cols(), |i, j| M.get(i, j))
}

pub fn inverse
<T: Num + FromPrimitive,
SM1: BasicMatrix<T>,
SM2: BasicMatrix<T>,
M: BasicMatrix<T> + SubMatrix<T, SM1>,
R: BasicMatrix<T> + Create<T> + SubMatrix<T, SM2>>(M: &M) -> R {
    // This basically does the blockwise inversion algorithm on the
    // Wikipedia page [1]. It's not a very efficient implementation,
    // since it ends up doing an absurd number of copies. It also
    // recurs down to 1x1 matrices, and in reality it makes sense to
    // do the cutoff sooner.
    //
    // [1] http://en.wikipedia.org/wiki/Invertible_matrix#Blockwise_inversion

    assert!(M.num_rows() == M.num_cols());

    let N = M.num_rows();

    do tracefn(("inverse", N)) {

    if N == 1 {
        Create::<T>::create(1, 1, |i, j| num::from_int::<T>(1).unwrap() / M.get(i, j))
    }
    else {
        let N2 = N / 2;
        let N2a = N - N2;

        // the convert is needed to avoid infinite expansion of
        // monomorphic functions. The type voodoo that's currently
        // commented out in matrix.rs might be able to remove the need
        // for this.
        let A: R = convert(&SubMatrix(M, 0, 0, N2, N2));
        let B = SubMatrix(M, 0,  N2, N2,  N2a);
        let C = SubMatrix(M, N2, 0,  N2a, N2);
        let D = SubMatrix(M, N2, N2, N2a, N2a);

        let Ai: R = inverse(&A);

        // Compute (D - CAiB)i
        //error!("Schur");
        let mut t: R = mat_mul(&C, &Ai);
        //error!("1");
        t = mat_mul(&t, &B);
        //error!("2");
        t = mat_sub(&D, &t);

        // new D
        //error!("D");
        let Dn: R = inverse(&t);

        // Compute AiB
        //error!("AiB");
        let AiB: R = mat_mul(&Ai, &B);
        
        // Compute CAi
        //error!("CAi");
        let CAi: R = mat_mul(&C, &Ai);

        // new C
        //error!("C");
        let mut Cn: R = mat_mul(&Dn, &CAi);
        mat_x_inplace(&mut Cn, num::from_int(-1).unwrap());

        // new B
        //error!("B");
        let mut Bn: R = mat_mul(&AiB, &Dn);
        // save the multiplication by -1 until later, since we can
        // reuse this result for the new A.

        // new A
        //error!("A");
        let mut An: R = mat_mul(&Bn, &CAi);
        An = mat_add(&Ai, &An);

        mat_x_inplace(&mut Bn, num::from_int(-1).unwrap());

        // Stitch it all back together.
        let top: R = concat_cols(&An, &Bn);
        let bot: R = concat_cols(&Cn, &Dn);
        concat_rows(&top, &bot)
    }
    }
}

pub fn cholesky_blocked
<SM: BasicMatrix<f64>, 
M: BasicMatrix<f64> + SubMatrix<f64, SM>,
R: BasicMatrix<f64> + Create<f64> + SubMatrix<f64, SM>>
(M: &M) -> R {
    /*
    A recursive blocked Cholesky factorization.

              +-----+
    Given M = | A B |
              | C D |
              +-----+

    The Cholesky decomposition should be... (c indicates Cholesky
    Factorization, i indicates inverse, t indicates transpose)
    
    +--------------------------+
    |  Ac            0         |
    |                          |
    | CActi  D - CActi(CActi)t |
    +--------------------------+

    I derived this with some help from
    http://www.netlib.org/utk/papers/factor/node9.html

    */

    assert!(M.num_rows() == M.num_cols());
    let N = M.num_rows();

    do tracefn(("cholesky_blocked", N)) {

    static BLOCK_SIZE: uint = 1;

    if N <= BLOCK_SIZE {
        let mut M = convert(M);
        cholesky_seq_inplace::<R>(&mut M);
        M
    }
    else {
        let N2 = N / 2;
        let N2a = N - N2;

        let A: R = convert(&SubMatrix(M, 0, 0, N2, N2));
        //let B = SubMatrix(M, 0,  N2, N2,  N2a);
        let C = SubMatrix(M, N2, 0,  N2a, N2);
        let D = SubMatrix(M, N2, N2, N2a, N2a);

        let Ac: R = cholesky_blocked(&A);

        let Aci: R = inverse::<f64, SubMatrix_t<f64,TransposeMatrix<f64,R>>, SM, TransposeMatrix<f64, R>, R>(&TransposeMatrix(&Ac));

        let CAci: R = mat_mul(&C, &Aci);

        let mut Dn: R = mat_mul::<f64, R, TransposeMatrix<f64, R>, R>(&CAci, &TransposeMatrix(&CAci));
        Dn = mat_sub(&D, &Dn);
        Dn = cholesky_blocked(&Dn);

        let Z = zero_matrix::<f64, R>(N2, N2a);

        let top: R = concat_cols(&Ac, &Z);
        let bot: R = concat_cols(&CAci, &Dn);

        concat_rows(&top, &bot)
    }
    }
}
