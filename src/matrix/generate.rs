use std::num::from_int;
use std::rand;
use std::rand::Rand;

use matrix::{BasicMatrix, Create, Matrix, Ring};

pub fn identity<T: Ring, M: BasicMatrix<T> + Create<T>>(N: uint)
    -> M
{
    do Create::<T>::create(N, N) |i, j| {
        if i == j {
            Ring::one()
        }
        else {
            Ring::zero()
        }
    }
}

// Generate a random square lower triangular matrix with unit diagonal.
pub fn rand_L1<T: Rand + FromPrimitive, M: BasicMatrix<T> + Create<T>>(N: uint)
    -> M {
    let mut r = rand::rng();

    do Create::create(N, N) |i, j| {
        if i == j {
            from_int::<T>(1).unwrap()
        }
        else if i > j {
            Rand::rand(&mut r)
        }
        else {
            from_int(0).unwrap()
        }
    }
}

pub fn zero_matrix<T: Ring, M: BasicMatrix<T> + Create<T>>
    (n: uint, m: uint) -> M
{
    Create::<T>::create(n, m, |_i, _j| Ring::zero())
}
