extern crate sync;
extern crate time;

extern crate SciRust = "SciRust#0.3pre";

use sync::Arc;

use SciRust::matrix;

use time::precise_time_s;

use SciRust::matrix::{BasicMatrix, Matrix, TransposeMatrix, transpose};
use SciRust::matrix::generate::{rand_L1};
use SciRust::matrix::algorithms::{mat_mul, cholesky_seq_inplace,
                        inverse, cholesky_blocked, mat_mul_blocked,
                        convert};
use SciRust::matrix::par;

// We'll settle for this for now.
//#[path="matrix/matrix.rs"]
//mod matrix;

type M = Matrix<f64>;

fn benchmark(N: uint) {
    println!("Benchmarking {:?} x {:?} matrices.", N, N);

    let L = rand_L1(N);
    let Lt = transpose(&L);

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
    let _Ab: M = mat_mul_blocked(&L, &Lt);
    let stop = precise_time_s();
    
    println!("Matrix Multiply (blocked): {:?}s", stop - start);

    let Ls: M = convert(&L);
    let Ls = Arc::new(Ls);
    let Lts: M = convert(&Lt);
    let Lts = Arc::new(Lts);

    //let start = precise_time_s();
    //println!("{:?}", ((*Ls).num_rows(), (*Lts).num_cols()));
    //let _Ap: M = par::mat_mul(&Ls, &Lts);
    //let stop = precise_time_s();
    //println!("Matrix Multiply (parallel): {:?}s", stop - start);

    let start = precise_time_s();
    let _Ai: M = inverse(&A);
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
    let _Ac: M = cholesky_blocked(&A);
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
