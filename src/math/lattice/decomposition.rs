use crate::math::component::big_matrix::BigMatrix;
use crate::math::component::big_vector::OwnedBigVector;
use malachite::base::num::arithmetic::traits::Abs;
use malachite::base::num::basic::traits::Zero;
use malachite::rational::Rational;

pub fn inverse(matrix: &BigMatrix) -> BigMatrix {
    let mut m = matrix.clone();
    let size = m.row_count;
    let mut p = OwnedBigVector::new(size);
    let mut inv = BigMatrix::identity(size);

    // Decomposition
    for i in 0..size {
        let mut pivot = size;
        let mut biggest_number = Rational::ZERO;

        for row in i..size {
            let d = m.get(row, i).abs();

            if d > biggest_number {
                biggest_number = d;
                pivot = row;
            }
        }

        if pivot == size {
            unreachable!("Matrix is singular");
        }

        p.set(i, Rational::from(pivot));
        inv.swap_rows(i, pivot);

        if pivot != i {
            m.swap_rows(i, pivot);
            // swaps += 1;
        }

        for row in i + 1..size {
            m.set(row, i, m.get(row, i) / m.get(i, i));
        }

        for row in i + 1..size {
            for col in i + 1..size {
                m.set(row, col, m.get(row, col) - (m.get(row, i) * m.get(i, col)));
            }
        }
    }

    // Inverse
    for dcol in 0..size {
        for row in 0..size {
            for col in 0..row {
                inv.set(
                    row,
                    dcol,
                    inv.get(row, dcol) - (m.get(row, col) * inv.get(col, dcol)),
                );
            }
        }
    }

    for dcol in 0..size {
        for row in (0..size).rev() {
            for col in (row + 1..size).rev() {
                inv.set(
                    row,
                    dcol,
                    inv.get(row, dcol) - (m.get(row, col) * inv.get(col, dcol)),
                );
            }

            inv.set(row, dcol, inv.get(row, dcol) / m.get(row, row));
        }
    }

    inv
}
