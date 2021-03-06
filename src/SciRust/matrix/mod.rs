// Linear Algrebra library for Rust
#![allow(non_snake_case_functions)]

use std::mem;
use std::num;
use std::ops;

pub mod algorithms;
pub mod generate;
pub mod util;
pub mod par;

// Sort of corresponds to the abstract algebra notion of a ring.
pub trait Ring : ops::Add<Self, Self> + ops::Mul<Self, Self> + FromPrimitive{
    fn one() -> Self;
    fn zero() -> Self;
}

impl<T: Num + FromPrimitive> Ring for T {
    fn one() -> T { num::from_int(1).unwrap() }
    fn zero() -> T { num::from_int(0).unwrap() }
}

pub trait BasicMatrix<T> {
    /// Returns an element in the matrix.
    fn get(&self, uint, uint) -> T;
    /// Set an element in the matrix.
    fn set(&mut self, uint, uint, T);

    /// Returns the number of rows in the matrix.
    fn num_rows(&self) -> uint;
    /// Returns the number of columns in the matrix.
    fn num_cols(&self) -> uint;

    #[lang="index"]
    fn index(&self, ix: &(uint, uint)) -> T {
        let &(i, j) = ix;
        self.get(i, j)
    }
}

pub trait Create<T> {
    fn create(uint, uint, |uint, uint| -> T) -> Self;
}

pub trait Vector<T> : ops::Index<uint, T> {
    fn len(&self) -> uint;
    fn get(&self, uint) -> T;
    fn set(&mut self, uint, T);
}

// Row and Column Vectors (Views into existing matrices)
pub struct RowVector<'r, T, M> {
    i: uint,
    base: &'r M
}

impl<'r, T, M: BasicMatrix<T>> Vector<T> for RowVector<'r, T, M> {
    fn len(&self) -> uint { self.base.num_cols() }
    #[inline(always)]
    fn get(&self, j: uint) -> T { self.base.get(self.i, j) }
    #[inline(always)]
    fn set(&mut self, j: uint, x: T) { self.base.set(self.i, j, x) }
}

impl<'r, T, M: BasicMatrix<T>> ops::Index<uint, T> for RowVector<'r, T, M> {
    #[inline(always)]
    fn index(&self, i: &uint) -> T { self.get(*i) }
}

pub struct ColumnVector<'r, T, M> {
    j: uint,
    base: &'r M
}

impl<'r, T, M: BasicMatrix<T>> Vector<T> for ColumnVector<'r, T, M> {
    fn len(&self) -> uint { self.base.num_rows() }
    #[inline(always)]
    fn get(&self, i: uint) -> T { self.base.get(i, self.j) }
    #[inline(always)]
    fn set(&mut self, i: uint, x: T) { self.base.set(i, self.j, x) }
}

impl<'r, T, M: BasicMatrix<T>> ops::Index<uint, T> for ColumnVector<'r, T, M> {
    #[inline(always)]
    fn index(&self, i: &uint) -> T { self.get(*i) }
}

impl<'r, T, M: BasicMatrix<T>> BasicMatrix<T> for &'r M {
    fn get(&self, i: uint, j: uint) -> T {
        self.get(i, j)
    }
    fn set(&mut self, i: uint, j: uint, x: T) {
        self.set(i, j, x)
    }

    fn num_rows(&self) -> uint {
        self.num_rows()
    }
    fn num_cols(&self) -> uint {
        self.num_rows()
    }
}

pub fn row<'a, T, M: BasicMatrix<T>>(m: &'a M, i: uint)
    -> RowVector<'a, T, M>
{
    assert!(i < m.num_rows());
    RowVector { i: i, base: m }
}

pub fn col<'a, T, M: BasicMatrix<T>>(m: &'a M, j: uint)
    -> ColumnVector<'a, T, M>
{
    assert!(j < m.num_cols());
    ColumnVector { j: j, base: m }
}

// A matrix in Row-Major Order
#[deriving(Clone)]
pub struct Matrix<T> {
    rows: uint,
    cols: uint,

    data: Vec<T>
}

impl<T: Clone> BasicMatrix<T> for Matrix<T> {
    #[inline(always)]
    fn get(&self, i: uint, j: uint) -> T {
        if i < self.num_rows() && j < self.num_cols() {
            self.data.as_slice()[i * self.num_cols() + j].clone()
        }
        else {
            fail!(format!("Index out of bounds. Index: {:?}, Dimension: {:?}",
                       (i, j),
                       (self.num_rows(), self.num_cols())))
        }
    }

    #[inline(always)]
    fn set(&mut self, i: uint, j: uint, x: T) {
        if i < self.num_rows() && j < self.num_cols() {
            let k = i * self.num_cols() + j;
            self.data.as_mut_slice()[k] = x
        }
        else {
            fail!(format!("Index out of bounds. Index: {:?}, Dimension: {:?}",
                       (i, j),
                       (self.num_rows(), self.num_cols())))
        }
    }

    fn num_rows(&self) -> uint { self.rows }
    fn num_cols(&self) -> uint { self.cols }
}

impl<T: Clone> Create<T> for Matrix<T> {
    fn create(i: uint, j: uint, init: |uint, uint| -> T)
        -> Matrix<T>
    {
        Matrix {
            rows: i,
            cols: j,
            data: Vec::from_fn(i * j, |k| {
                let i = k / j;
                let j = k % j;
                init(i, j)
            })
        }
    }
}

pub struct SubMatrix<'r, T, M> {
    i: uint, j: uint,
    rows: uint, cols: uint,
    base: &'r M
}

pub fn SubMatrix<'a, T, M: BasicMatrix<T>>(m: &'a M,
                                         i: uint,
                                         j: uint,
                                         rows: uint,
                                         cols: uint)
    -> SubMatrix<'a, T, M>
{
    assert!(rows > 0);
    assert!(cols > 0);
    assert!(i < m.num_rows());
    assert!(j < m.num_cols());
    assert!(i + rows <= m.num_rows());
    assert!(j + cols <= m.num_cols());
    SubMatrix {
        i: i, j: j, rows: rows, cols: cols, base: m
    }        
}

impl<'r, T, M: BasicMatrix<T>> BasicMatrix<T> for SubMatrix<'r, T, M> {
    fn num_rows(&self) -> uint { self.rows }
    fn num_cols(&self) -> uint { self.cols }

    #[inline(always)]
    fn get(&self, i: uint, j: uint) -> T {
        if i < self.rows && j < self.cols {
            self.base.get(i + self.i, j + self.j)
        }
        else {
            fail!("SubMatrix index out of bounds.")
        }
    }

    #[inline(always)]
    fn set(&mut self, i: uint, j: uint, x: T) {
        if i < self.rows && j < self.cols {
            self.base.set(i + self.i, j + self.j, x)
        }
        else {
            fail!("SubMatrix index out of bounds.")
        }
    }
}

pub struct TransposeMatrix<'r, T, M>(&'r M);

impl<'r, T, M> Clone for TransposeMatrix<'r, T, M> {
    fn clone(&self) -> TransposeMatrix<'r, T, M> {
        match self {
            &TransposeMatrix(m) => TransposeMatrix(m)
        }
    }
}

pub fn transpose<'r, T, M>(m: &'r M) -> TransposeMatrix<'r, T, M> {
    TransposeMatrix(m)
}

impl<'r, T, M> TransposeMatrix<'r, T, M> {
    fn get_ref(&self) -> &'r M {
        match *self {
            TransposeMatrix(m) => m
        }
    }

    fn get_mut(&mut self) -> &'r mut M {
        match self {
            &TransposeMatrix(m) => unsafe { mem::transmute(m) }
        }
    }
}

impl<'r, T, M: BasicMatrix<T>> BasicMatrix<T> for TransposeMatrix<'r, T, M> {
    fn num_rows(&self) -> uint { 
        self.get_ref().num_cols()
    }
    fn num_cols(&self) -> uint {
        self.get_ref().num_rows()
    }

    #[inline(always)]
    fn get(&self, i: uint, j: uint) -> T {
        self.get_ref().get(j, i)
    }

    #[inline(always)]
    fn set(&mut self, i: uint, j: uint, x: T) {
        self.get_mut().set(j, i, x)
    }
}
