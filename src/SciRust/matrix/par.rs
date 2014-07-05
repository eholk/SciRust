#![allow(uppercase_variables)]

use std::sync::Arc;
use std::sync::Future;

use matrix::{BasicMatrix, Create, SubMatrix};
use matrix::algorithms;
use matrix::algorithms::{concat_cols, concat_rows, mat_add_inplace};

type SubCoords = ((uint, uint), (uint, uint));

fn subdivide_coords(x: SubCoords) -> (SubCoords, SubCoords, SubCoords, SubCoords) {
    let ((i, j), (n, m)) = x;
    let n2 = n / 2;
    let m2 = m / 2;
    let n2a = n - n2;
    let m2a = m - m2;

    (((i, j), (n2, m2)),
     ((i + n2, j), (n2a, m2)),
     ((i, j + m2), (n2, m2a)),
     ((i + n2, j + m2), (n2a, m2a)))
}

pub fn mat_mul<T: Num + FromPrimitive, LHS: BasicMatrix<T> + Send + Clone, RHS: BasicMatrix<T> + Send + Clone, Res: BasicMatrix<T> + Create<T> + Send>
(lhs: &LHS, rhs: &RHS) -> Res {
    assert!(lhs.num_cols() == rhs.num_rows());

    sub_mul(lhs.clone(), ((0, 0), (lhs.num_rows(), lhs.num_cols())),
            rhs.clone(), ((0, 0), (rhs.num_rows(), rhs.num_cols())))
}

fn sub_mul<T: Num + FromPrimitive, LHS: BasicMatrix<T> + Send + Clone, RHS: BasicMatrix<T> + Send + Clone, Res: BasicMatrix<T> + Create<T> + Send>
(lhs: LHS, lc: SubCoords, rhs: RHS, rc: SubCoords) -> Res
{
    static BLOCK_SIZE: uint = 1 << 13;

    let ((li, lj), (ln, lm)) = lc;
    let ((ri, rj), (rn, rm)) = rc;
    if ln * lm <= BLOCK_SIZE || rn * rm <= BLOCK_SIZE {
        algorithms::mat_mul(&SubMatrix(&lhs, li, lj, ln, lm),
                            &SubMatrix(&rhs, ri, rj, rn, rm))
    }
        else {
        let (A, B, C, D) = subdivide_coords(lc);
        let (E, F, G, H) = subdivide_coords(rc);
        
        let lhsc = lhs.clone();
        let rhsc = rhs.clone();
        let AE: Future<Res> 
            = Future::spawn(proc() {
                sub_mul(lhsc.clone(), A, rhsc.clone(), E)
            });
        let lhsc = lhs.clone();
        let rhsc = rhs.clone();
        let BG: Future<Res> 
            = Future::spawn(proc() {
                sub_mul(lhsc.clone(), B, rhsc.clone(), G)
            });
        let lhsc = lhs.clone();
        let rhsc = rhs.clone();
        let AF: Future<Res> 
            = Future::spawn(proc() {
                sub_mul(lhsc.clone(), A, rhsc.clone(), F)
            });
        let lhsc = lhs.clone();
        let rhsc = rhs.clone();
        let BH: Future<Res> 
            = Future::spawn(proc() {
                sub_mul(lhsc.clone(), B, rhsc.clone(), H)
            });

        let lhsc = lhs.clone();
        let rhsc = rhs.clone();
        let CE: Future<Res> 
            = Future::spawn(proc() {
                sub_mul(lhsc.clone(), C, rhsc.clone(), E)
            });
        let lhsc = lhs.clone();
        let rhsc = rhs.clone();
        let DG: Future<Res> 
            = Future::spawn(proc() {
                sub_mul(lhsc.clone(), D, rhsc.clone(), G)
            });
        let lhsc = lhs.clone();
        let rhsc = rhs.clone();
        let CF: Future<Res> 
            = Future::spawn(proc() {
                sub_mul(lhsc.clone(), C, rhsc.clone(), F)
            });
        let DH: Future<Res> 
            = Future::spawn(proc() sub_mul(lhs.clone(), D, rhs.clone(), H));
        
        let mut An = AE.unwrap();
        mat_add_inplace(&mut An, &BG.unwrap());

        let mut Bn = AF.unwrap();
        mat_add_inplace(&mut Bn, &BH.unwrap());

        let mut Cn = CE.unwrap();
        mat_add_inplace(&mut Cn, &DG.unwrap());

        let mut Dn = CF.unwrap();
        mat_add_inplace(&mut Dn, &DH.unwrap());

        let top: Res = concat_cols(&An, &Bn);
        let bot: Res = concat_cols(&Cn, &Dn);
        concat_rows(&top, &bot)
    }
}

impl<T, M: BasicMatrix<T> + Send> BasicMatrix<T> for Arc<M> {
    fn get(&self, i: uint, j: uint) -> T {
        (*self).get(i, j)
    }

    fn set(&mut self, _i: uint, _j: uint, _x: T) {
        fail!("Attempting to mutate shared matrix.");
    }

    fn num_rows(&self) -> uint { (*self).num_rows() }
    fn num_cols(&self) -> uint { (*self).num_cols() }
}
