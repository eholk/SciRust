use num::Num;

fn identity<T: Copy Ring, M: BasicMatrix<T> Create<T, M>>(N: uint)
    -> M
{
    // Why does create have three type parameters?
    do create::<T, M, M>(N, N) |i, j| {
        if i == j {
            one::<T>()
        }
        else {
            zero()
        }
    }
}

// Generate a random square lower triangular matrix with unit diagonal.
fn rand_L1(N: uint) -> Matrix<float> {
    let r = rand::Rng();
    do create::<float, Matrix<float>, Matrix<float>>(N, N) |i, j| {
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

fn zero_matrix<T: Copy Ring, M: BasicMatrix<T> Create<T, M>>(n: uint, m: uint) -> M
{
    create::<T, M, M>(n, m, |_i, _j| zero::<T>())
}