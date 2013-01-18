use num::Num;

use matrix::{BasicMatrix, Create, Matrix, Ring};

pub fn identity<T: Copy Ring, M: BasicMatrix<T> Create<T, M>>(N: uint)
    -> M
{
    // Why does create have three type parameters?
    do Create::create::<T, M, M>(N, N) |i, j| {
        if i == j {
            Ring::one::<T>()
        }
        else {
            Ring::zero()
        }
    }
}

// Generate a random square lower triangular matrix with unit diagonal.
pub fn rand_L1(N: uint) -> Matrix<float> {
    let r = rand::Rng();
    do Create::create::<float, Matrix<float>, Matrix<float>>(N, N) |i, j| {
        if i == j {
            1.0
        }
        else if i > j {
            r.gen_float()
        }
        else {
            0f
        }
    }
}

pub fn zero_matrix<T: Copy Ring, M: BasicMatrix<T> Create<T, M>>(n: uint, m: uint) -> M
{
    Create::create::<T, M, M>(n, m, |_i, _j| Ring::zero::<T>())
}
