extern mod extra;

// This version makes ICE
//extern mod SciRust;
//use SciRust::matrix;

use extra::time::precise_time_s;

use matrix::{Matrix, TransposeMatrix, Create};
use matrix::generate::{identity, rand_L1};
use matrix::algorithms::{dot, mat_mul, transpose, cholesky_seq_inplace,
                        inverse, cholesky_blocked, mat_mul_blocked,
                        convert};
use matrix::util::to_str;

// We'll settle for this for now.
#[path="matrix/matrix.rs"]
mod matrix;

type M = Matrix<f64>;

fn benchmark(N: uint) {
    println!("Benchmarking {:?} x {:?} matrices.", N, N);

    let L = rand_L1(N);
    let Lt = TransposeMatrix(&L);

    let start = precise_time_s();
    let A: M = mat_mul::<f64, M, TransposeMatrix<f64, M>, M>(&L, &Lt);
    let stop = precise_time_s();

    println!("Matrix Multiply: {:?}s", stop - start);

    //let start = precise_time_s();
    //let Ap: M = par::mat_mul(&L, &Lt);
    //let stop = precise_time_s();
    //
    //println!("Matrix Multiply (parallel): {:?}s", stop - start);

    // TODO: make sure A and Ap agree.

    let start = precise_time_s();
    let Ab: M = mat_mul_blocked(&L, &Lt);
    let stop = precise_time_s();

    println!("Matrix Multiply (blocked): {:?}s", stop - start);

    let start = precise_time_s();
    let Ai: M = inverse(&A);
    let stop = precise_time_s();

    println!("Matrix Inverse: {:?}s", stop - start);

    //let start = precise_time_s();
    //let Ai: M = par::inverse(&A);
    //let stop = precise_time_s();
    //
    //println!("Matrix Inverse (parallel): {:?}s", stop - start);

    let mut A2 = convert(&A);
    let start = precise_time_s();
    cholesky_seq_inplace::<M>(&mut A2);
    let stop = precise_time_s();

    println!("Cholesky (sequential): {:?}s", stop - start);

    let start = precise_time_s();
    let Ac: M = cholesky_blocked(&A);
    let stop = precise_time_s();

    println!("Cholesky (blocked): {:?}s", stop - start);   

    //let start = precise_time_s();
    //let Ac: M = par::cholesky_blocked(&A);
    //let stop = precise_time_s();
    //
    //println!("Cholesky (parallel): {:?}s", stop - start);   
}

fn main() {
    //benchmark(16);
    benchmark(1200);
}
