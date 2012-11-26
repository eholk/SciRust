use matrix::generate::zero_matrix;

pub mod par;

pub fn dot<T: Copy Num, L: Vector<T>, R: Vector<T>>(lhs: &L, rhs: &R) -> T {
    if lhs.len() != rhs.len() {
        fail ~"Invalid vector lengths."
    }

    //error!("%? ### %?", lhs, rhs);

    let mut acc : T = num::from_int(0);
    for uint::range(0, lhs.len()) |i| {
        acc += lhs[i] * rhs[i]
    }

    acc
}

pub fn mat_mul<T: Copy Num, LHS: BasicMatrix<T>, RHS: BasicMatrix<T>, Res: BasicMatrix<T> Create<T, Res>> (lhs: &LHS, rhs: &RHS) -> Res
{
    if lhs.num_cols() != rhs.num_rows() {
        fail fmt!("Incompatible matrix sizes. LHS: %?, RHS: %?",
                  (lhs.num_rows(), lhs.num_cols()),
                  (rhs.num_rows(), rhs.num_cols()))
    }

    //error!("Multiplying %? by %? -> %?",
    //       (lhs.num_rows(), lhs.num_cols()),
    //       (rhs.num_rows(), rhs.num_cols()),
    //       (lhs.num_rows(), rhs.num_cols()));

    do create::<T, Res, Res>(lhs.num_rows(), rhs.num_cols()) |i, j| {
        dot(&row(lhs, i), &col(rhs, j))
    }
}

pub fn mat_add<T: Copy Num, LHS: BasicMatrix<T>, RHS: BasicMatrix<T>, Res: BasicMatrix<T> Create<T, Res>> (lhs: &LHS, rhs: &RHS) -> Res
{
    if lhs.num_cols() != rhs.num_cols() || lhs.num_rows() != rhs.num_rows() {
        fail ~"Incompatible matrix sizes"
    }

    do create::<T, Res, Res>(lhs.num_rows(), rhs.num_cols()) |i, j| {
        lhs.get(i, j) + rhs.get(i, j)
    }
}

pub fn mat_sub<T: Copy Num, LHS: BasicMatrix<T>, RHS: BasicMatrix<T>, Res: BasicMatrix<T> Create<T, Res>> (lhs: &LHS, rhs: &RHS) -> Res
{
    if lhs.num_cols() != rhs.num_cols() || lhs.num_rows() != rhs.num_rows() {
        fail fmt!("Incompatible matrix sizes. LHS: %?, RHS: %?",
                  (lhs.num_rows(), lhs.num_cols()),
                  (rhs.num_rows(), rhs.num_cols()))
    }

    do create::<T, Res, Res>(lhs.num_rows(), rhs.num_cols()) |i, j| {
        lhs.get(i, j) - rhs.get(i, j)
    }
}

pub fn mat_x_inplace<T: Copy Num, M: BasicMatrix<T>>(A: &M, x: T) {
    for_each(A, |_i, _j, y| x * y);
}

pub fn transpose<T: Copy, M: BasicMatrix<T>, R: BasicMatrix<T> Create<T, R>>(m: &M) -> R {
    do create::<T, R, R>(m.num_cols(), m.num_rows()) |i, j| {
        m.get(j, i)
    }
}

pub fn cholesky_seq_inplace_raw<T: Copy Num Sqrt, M: BasicMatrix<float>>(A: &M, start: uint) {
    assert A.num_rows() == A.num_cols();
    let N = A.num_rows();
    for uint::range(start, N) |k| {
        let Akk = A.get(k, k);

        A.set(k, k, Akk.sqrt());

        for uint::range(k + 1, N) |i| {
            let Aik = A.get(i, k);
            A.set(i, k, Aik / Akk);
        }

        for uint::range(k + 1, N) |i| {
            let Aik = A.get(i, k);
            for uint::range(k + 1, i + 1) |j| {
                let Ajk = A.get(j, k);
                let Aij = A.get(i, j);
                A.set(i, j, Aij - Aik * Ajk);
            }
        }
    }
}

pub fn for_each<T: Copy, M: BasicMatrix<T>>(A: &M, f: fn(uint, uint, T) -> T) {
    for uint::range(0, A.num_rows()) |i| {
        for uint::range(0, A.num_cols()) |j| {
            A.set(i, j, f(i, j, A.get(i, j)))
        }
    }
}

pub fn cholesky_seq_inplace<T: Copy Num Sqrt, M: BasicMatrix<float>>(A: &M) {
    cholesky_seq_inplace_start::<T, M>(A, 0);
}

pub fn cholesky_seq_inplace_start<T: Copy Num Sqrt, M: BasicMatrix<float>>(A: &M,
                                                                 start: uint) {
    cholesky_seq_inplace_raw::<T, M>(A, start);

    for uint::range(start, A.num_rows()) |i| {
        for uint::range(i + 1, A.num_cols()) |j| {
            A.set(i, j, num::from_int(0))
        }
    }    
}


pub fn concat_rows<T: Copy, LHS: BasicMatrix<T>, RHS: BasicMatrix<T>, R: BasicMatrix<T> Create<T, R>>(A: &LHS, B: &RHS) -> R {
    assert A.num_cols() == B.num_cols();

    let N = A.num_rows();

    do create::<T, R, R>(N + B.num_rows(), A.num_cols()) |i, j| {
        if i < N {
            A.get(i, j)
        }
        else {
            B.get(i - N, j)
        }
    }
}

pub fn concat_cols<T: Copy, LHS: BasicMatrix<T>, RHS: BasicMatrix<T>, R: BasicMatrix<T> Create<T, R>>(A: &LHS, B: &RHS) -> R {
    //error!("concat: %?, %?",
    //       (A.num_rows(), A.num_cols()),
    //       (B.num_rows(), B.num_cols()));
           
    assert A.num_rows() == B.num_rows();

    let N = A.num_cols();

    do create::<T, R, R>(B.num_rows(), N + B.num_cols()) |i, j| {
        if j < N {
            A.get(i, j)
        }
        else {
            B.get(i, j - N)
        }
    }
}

pub fn convert<T: Copy, M: BasicMatrix<T>, R: BasicMatrix<T> Create<T, R>>(M: &M) -> R {
    create::<T, R, R>(M.num_rows(), M.num_cols(), |i, j| M.get(i, j))
}

pub fn inverse<T: Copy Num, M: BasicMatrix<T>, R: BasicMatrix<T> Create<T, R>>(M: &M) -> R {
    // This basically does the blockwise inversion algorithm on the
    // Wikipedia page [1]. It's not a very efficient implementation,
    // since it ends up doing an absurd number of copies. It also
    // recurs down to 1x1 matrices, and in reality it makes sense to
    // do the cutoff sooner.
    //
    // [1] http://en.wikipedia.org/wiki/Invertible_matrix#Blockwise_inversion

    assert M.num_rows() == M.num_cols();

    let N = M.num_rows();

    if N == 1 {
        create::<T, R, R>(1, 1, |i, j| num::from_int::<T>(1) / M.get(i, j))
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
        let Cn: R = mat_mul(&Dn, &CAi);
        mat_x_inplace(&Cn, num::from_int(-1));

        // new B
        //error!("B");
        let Bn: R = mat_mul(&AiB, &Dn);
        // save the multiplication by -1 until later, since we can
        // reuse this result for the new A.

        // new A
        //error!("A");
        let mut An: R = mat_mul(&Bn, &CAi);
        An = mat_add(&Ai, &An);

        mat_x_inplace(&Bn, num::from_int(-1));

        // Stitch it all back together.
        let top: R = concat_cols(&An, &Bn);
        let bot: R = concat_cols(&Cn, &Dn);
        concat_rows(&top, &bot)
    }
}

pub fn cholesky_blocked<M: BasicMatrix<float>, R: BasicMatrix<float> Create<float, R>>(M: &M) -> R {
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

    assert M.num_rows() == M.num_cols();
    let N = M.num_rows();

    const BLOCK_SIZE: uint = 1;

    if N <= BLOCK_SIZE {
        let M = convert(M);
        cholesky_seq_inplace::<float, R>(&M);
        move M
    }
    else {
        let N2 = N / 2;
        let N2a = N - N2;

        let A: R = convert(&SubMatrix(M, 0, 0, N2, N2));
        //let B = SubMatrix(M, 0,  N2, N2,  N2a);
        let C = SubMatrix(M, N2, 0,  N2a, N2);
        let D = SubMatrix(M, N2, N2, N2a, N2a);

        let Ac: R = cholesky_blocked(&A);

        let Aci: R = inverse(&TransposeMatrix(&Ac));

        let CAci: R = mat_mul(&C, &Aci);

        let mut Dn: R = mat_mul(&CAci, &TransposeMatrix(&CAci));
        Dn = mat_sub(&D, &Dn);
        Dn = cholesky_blocked(&Dn);

        let Z = zero_matrix::<float, R>(N2, N2a);

        let top: R = concat_cols(&Ac, &Z);
        let bot: R = concat_cols(&CAci, &Dn);

        concat_rows(&top, &bot)
    }
}
