// Linear Algrebra library for Rust

pub mod algorithms;
pub mod generate;
pub mod util;

use num::Num;
use to_str::ToStr;

// Sort of corresponds to the abstract algebra notion of a ring.
//
// FIXME: This should also include Copy, Add and Mul, but it doesn't
// due to bugs in Rust.
pub trait Ring : {
    static fn one() -> self;
    static fn zero() -> self;
}

pub impl<T: num::Num> T : Ring {
    static fn one() -> T { Num::from_int(1) }
    static fn zero() -> T { Num::from_int(0) }
}

pub trait BasicMatrix<T: Copy> {
    pure fn get(uint, uint) -> T;
    fn set(uint, uint, T);

    pure fn num_rows() -> uint;
    pure fn num_cols() -> uint;
}

pub trait Create<T: Copy, M: BasicMatrix<T>> {
    static fn create(uint, uint, fn(uint, uint) -> T) -> M;
}

pub trait Vector<T: Copy> : ops::Index<uint, T> {
    pure fn len(&self) -> uint;
    pure fn get(uint) -> T;
    fn set(uint, T);
}

// It's nice to be able to use arbitrary vector slices as Vectors.
pub impl<T: Copy> &[mut T] : Vector<T> {
    pure fn len(&self) -> uint { vec::len(*self) }
    #[inline(always)]
    pure fn get(i: uint) -> T { self[i] }
    #[inline(always)]
    fn set(i: uint, x: T) { self[i] = x }
}

pub impl<T: Copy> &[mut T] : ops::Index<uint, T> {
    #[inline(always)]
    pure fn index(&self, i: uint) -> T { self[i] }
}


pub impl<T: Copy> ~[mut T] : Vector<T> {
    pure fn len(&self) -> uint { vec::len(*self) }
    #[inline(always)]
    pure fn get(i: uint) -> T { self[i] }
    #[inline(always)]
    fn set(i: uint, x: T) { self[i] = x }
}

pub impl<T: Copy> ~[mut T] : ops::Index<uint, T> {
    #[inline(always)]
    pure fn index(&self, i: uint) -> T { self[i] }
}

// Row and Column Vectors (Views into existing matrices)
struct RowVector<T: Copy, M: BasicMatrix<T>> {
    i: uint,
    base: &M
}

impl<T: Copy, M: BasicMatrix<T>> RowVector<T, M> : Vector<T> {
    pure fn len(&self) -> uint { self.base.num_cols() }
    #[inline(always)]
    pure fn get(j: uint) -> T { self.base.get(self.i, j) }
    #[inline(always)]
    fn set(j: uint, x: T) { self.base.set(self.i, j, x) }
}

impl<T: Copy, M: BasicMatrix<T>> RowVector<T, M> : ops::Index<uint, T> {
    #[inline(always)]
    pure fn index(&self, i: uint) -> T { self.get(i) }
}

struct ColumnVector<T: Copy, M: BasicMatrix<T>> {
    j: uint,
    base: &M
}

impl<T: Copy, M: BasicMatrix<T>> ColumnVector<T, M> : Vector<T> {
    pure fn len(&self) -> uint { self.base.num_rows() }
    #[inline(always)]
    pure fn get(i: uint) -> T { self.base.get(i, self.j) }
    #[inline(always)]
    fn set(i: uint, x: T) { self.base.set(i, self.j, x) }
}

impl<T: Copy, M: BasicMatrix<T>> ColumnVector<T, M> : ops::Index<uint, T> {
    #[inline(always)]
    pure fn index(&self, i: uint) -> T { self.get(i) }
}

pub pure fn row<T: Copy, M: BasicMatrix<T>>(m: &a/M, i: uint)
    -> RowVector/&a<T, M>
{
    RowVector { i: i, base: m }
}

pub pure fn col<T: Copy, M: BasicMatrix<T>>(m: &a/M, j: uint)
    -> ColumnVector/&a<T, M>
{
    ColumnVector { j: j, base: m }
}

// A matrix in Row-Major Order
pub struct Matrix/&<T: Copy> {
    rows: uint,
    cols: uint,

    data: ~[mut T]
}

pub impl<T: Copy> Matrix<T> : BasicMatrix<T> {
    #[inline(always)]
    pure fn get(i: uint, j: uint) -> T {
        if i < self.num_rows() && j < self.num_cols() {
            self.data[i * self.num_cols() + j]
        }
        else {
            fail fmt!("Index out of bounds. Index: %?, Dimension: %?",
                      (i, j),
                      (self.num_rows(), self.num_cols()))
        }
    }

    #[inline(always)]
    fn set(i: uint, j: uint, x: T) {
        if i < self.num_rows() && j < self.num_cols() {
            self.data[i * self.num_cols() + j] = x
        }
        else {
            fail fmt!("Index out of bounds. Index: %?, Dimension: %?",
                      (i, j),
                      (self.num_rows(), self.num_cols()))
        }
    }

    pure fn num_rows() -> uint { self.rows }
    pure fn num_cols() -> uint { self.cols }
}

pub impl<T: Copy> Matrix<T> : Create<T, Matrix<T>> {
    static fn create(i: uint, j: uint, init: fn(uint, uint) -> T)
        -> Matrix<T>
    {
        Matrix {
            rows: i,
            cols: j,
            data: vec::to_mut(do vec::from_fn(i * j) |k| {
                let i = k / j;
                let j = k % j;
                init(i, j)
            })
        }
    }
}

struct SubMatrix<T: Copy, M: BasicMatrix<T>> {
    i: uint, j: uint,
    rows: uint, cols: uint,
    base: &M
}

/* // This is some type voodoo that might be nice later.

trait RefineSubMatrix<T: Copy, M: BasicMatrix<T>> {
    fn refine_submatrix(m: &self,
                        i: uint, j: uint,
                        rows: uint, cols: uint) -> SubMatrix<T, M>;
}

impl<T: Copy> Matrix<T>: RefineSubMatrix<T, Matrix<T>> {
    fn refine_submatrix(m: &a/Matrix<T>,
                        i: uint, j: uint,
                        rows: uint, cols: uint)
        -> SubMatrix/&a<T, Matrix<T>>
    {
        SubMatrix {
            i: i, j: j, rows: rows, cols: cols, base: m
        }        
    }
}

impl<T: Copy, M: BasicMatrix<T>> SubMatrix<T, M>: RefineSubMatrix<T, M> {
    fn refine_submatrix(m: &a/SubMatrix<T, M>,
                        i: uint, j: uint,
                        rows: uint, cols: uint)
        -> SubMatrix/&a<T, M>
    {
        SubMatrix {
            i: i, j: j, rows: rows, cols: cols, base: m.base
        }        
    }    
}

fn SubMatrix<T: Copy, M: BasicMatrix<T> RefineSubMatrix<T, M>>(m: &a/M,
                                                               i: uint,
                                                               j: uint,
                                                               rows: uint,
                                                               cols: uint)
    -> SubMatrix/&a<T, M>
{
    m.refine_submatrix(m, i, j, rows, cols)
}
*/

pub fn SubMatrix<T: Copy, M: BasicMatrix<T>>(m: &a/M,
                                         i: uint,
                                         j: uint,
                                         rows: uint,
                                         cols: uint)
    -> SubMatrix/&a<T, M>
{
    SubMatrix {
        i: i, j: j, rows: rows, cols: cols, base: m
    }        
}

impl<T: Copy, M: BasicMatrix<T>> SubMatrix<T, M> : BasicMatrix<T> {
    pure fn num_rows() -> uint { self.rows }
    pure fn num_cols() -> uint { self.cols }

    #[inline(always)]
    pure fn get(i: uint, j: uint) -> T {
        if i < self.rows && j < self.cols {
            self.base.get(i + self.i, j + self.j)
        }
        else {
            fail ~"SubMatrix index out of bounds."
        }
    }

    #[inline(always)]
    fn set(i: uint, j: uint, x: T) {
        if i < self.rows && j < self.cols {
            self.base.set(i + self.i, j + self.j, x)
        }
        else {
            fail ~"SubMatrix index out of bounds."
        }
    }
}

pub trait Sqrt {
    pure fn sqrt() -> self;
}

impl float: Sqrt {
    pure fn sqrt() -> float {
        float::sqrt(self)
    }
}

enum TransposeMatrix<T: Copy, M: BasicMatrix<T>> = &M;

impl<T: Copy, M: BasicMatrix<T>> TransposeMatrix<T, M>: BasicMatrix<T> {
    pure fn num_rows() -> uint { (*self).num_cols() }
    pure fn num_cols() -> uint { (*self).num_rows() }

    #[inline(always)]
    pure fn get(i: uint, j: uint) -> T { (*self).get(j, i) }
    #[inline(always)]
    fn set(i: uint, j: uint, x: T) { (*self).set(j, i, x) }
}