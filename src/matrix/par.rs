extern mod std;

use num::Num;

use std::future::spawn;

use matrix::{BasicMatrix, SubMatrix, TransposeMatrix, Create, row, col};
use matrix::algorithms::{concat_rows, concat_cols, dot, convert, mat_sub,
                        mat_x_inplace, mat_add, cholesky_seq_inplace,
                        transpose};
use matrix::generate::{zero_matrix};

// A parallel matrix creator.
pub fn create<T: Copy Owned, M: Owned Copy BasicMatrix<T> Create<T, M>>(rows: uint, cols: uint, f: fn~(uint, uint) -> T) -> M {
    // by default, use 128x128 as a block size. We should use
    // autotuning to determine the best option.

    create_blocked(rows, cols, 128, move f)
}

fn min<T: cmp::Ord>(a: T, b: T) -> T {
    if a < b { move a } else { move b }
}

pub fn create_blocked<T: Copy Owned, M: Owned Copy BasicMatrix<T> Create<T, M>>(rows: uint, cols: uint, block_size: uint, f: fn~(uint, uint) -> T) -> M {
    
    let mut blocks = ~[];

    let mut i = 0;

    while i < rows {
        let mut row = ~[];
        let rows = min(block_size, rows - i);

        let mut j = 0;
        while j < cols {
            
            let cols = min(block_size, cols - j);

            let fu = do spawn |copy i, copy j, copy f| {
                //error!("block %?: %?", (i, j), (rows, cols));
                do Create::create::<T, M, M>(rows, cols) |ii, jj|
                {
                    f(i + ii, j + jj)
                }
            };

            row.push(move fu);
            j += block_size;
        }
        
        blocks.push(move row);
        i += block_size;
    }
    
    // read the first row
    let mut M: M = blocks[0][0].get();
    for uint::range(1, blocks[0].len()) |i| {
        M = concat_cols(&M, &blocks[0][i].get())
    }
    
    // read the rest of the rows
    for uint::range(1, blocks.len()) |j| {
        let mut r: M = blocks[j][0].get();
        for uint::range(1, blocks[j].len()) |i| {
            r = concat_cols(&r, &blocks[j][i].get())
        }
        
        M = concat_rows(&M, &r)
    }

    M
}

pub fn mat_mul<T: Owned Copy Num, LHS: BasicMatrix<T>, RHS: BasicMatrix<T>, Res: Owned Copy BasicMatrix<T> Create<T, Res>> (lhs: &LHS, rhs: &RHS) -> Res
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

    // In a perfect world, we'd use ARCs. For now, we'll use unsafety.
    unsafe {
        let num_rows = lhs.num_rows();
        let num_cols = rhs.num_cols();
        let lhs = ptr::addr_of(lhs);
        let rhs = ptr::addr_of(rhs);
        let M = do create::<T, Res>(num_rows, num_cols) |i, j| unsafe {
            let lhs: &LHS = &*lhs;
            let rhs: &RHS = &*rhs;
            dot(&row(lhs, i), &col(rhs, j))
        };

        //error!("Done");
        M
    }
}

pub fn inverse<T: Copy Owned Num, M: Owned BasicMatrix<T>, R: Copy Owned BasicMatrix<T> Create<T, R>>(M: &M) -> R {
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
        Create::create::<T, R, R>(1, 1, |i, j|
                                  Num::from_int::<T>(1) / M.get(i, j))
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
        mat_x_inplace(&Cn, Num::from_int(-1));

        // new B
        //error!("B");
        let Bn: R = mat_mul(&AiB, &Dn);
        // save the multiplication by -1 until later, since we can
        // reuse this result for the new A.

        // new A
        //error!("A");
        let mut An: R = mat_mul(&Bn, &CAi);
        An = mat_add(&Ai, &An);

        mat_x_inplace(&Bn, Num::from_int(-1));

        // Stitch it all back together.
        let top: R = concat_cols(&An, &Bn);
        let bot: R = concat_cols(&Cn, &Dn);
        concat_rows(&top, &bot)
    }
}

pub fn cholesky_blocked<M: Owned Copy BasicMatrix<float>, R: Owned Copy BasicMatrix<float> Create<float, R>>(M: &M) -> R {
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
        let C: R = convert(&SubMatrix(M, N2, 0,  N2a, N2));
        let D: R = convert(&SubMatrix(M, N2, N2, N2a, N2a));

        let Ac: R = cholesky_blocked(&A);

        let Act: R = transpose(&Ac);
        let Aci: R = inverse(&Act);

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
