use crate::math::component::big_vector::Dot;
use crate::math::component::big_vector::{MutViewBigVector, OwnedBigVector, ViewBigVector};
use crate::util::util::filled_vec;
use malachite::base::num::basic::traits::{One, Zero};
use malachite::rational::Rational;
use std::ops::{Add, Mul, Sub};

pub struct BigMatrix {
    pub numbers: Vec<Rational>,
    pub row_count: usize,
    pub column_count: usize,
}

impl BigMatrix {
    pub fn identity(size: usize) -> BigMatrix {
        let mut matrix = Self::new(size, size);

        for i in 0..size {
            matrix.set(i, i, Rational::ONE);
        }

        matrix
    }

    pub fn new(row_count: usize, column_count: usize) -> BigMatrix {
        BigMatrix {
            numbers: filled_vec(row_count * column_count, |_i| Rational::ZERO),
            row_count,
            column_count,
        }
    }

    pub fn new_provided<F: Fn(usize, usize) -> Rational>(
        row_count: usize,
        column_count: usize,
        generator: F,
    ) -> BigMatrix {
        let mut numbers: Vec<Rational> = Vec::with_capacity(row_count * column_count);

        for row in 0..row_count {
            for col in 0..column_count {
                numbers.push(generator(row, col));
            }
        }

        BigMatrix {
            numbers,
            row_count,
            column_count,
        }
    }

    pub fn is_square(&self) -> bool {
        self.row_count == self.column_count
    }

    #[inline]
    pub fn get(&self, row: usize, col: usize) -> &Rational {
        &self.numbers[col + self.column_count * row]
    }

    #[inline]
    pub fn set(&mut self, row: usize, col: usize, value: Rational) {
        self.numbers[col + self.column_count * row] = value;
    }

    pub fn get_row(&'_ self, row: usize) -> ViewBigVector<'_> {
        ViewBigVector::new(self.column_count, &self.numbers, row * self.column_count, 1)
    }

    pub fn get_row_mut(&'_ mut self, row: usize) -> MutViewBigVector<'_> {
        MutViewBigVector::new(
            self.column_count,
            &mut self.numbers,
            row * self.column_count,
            1,
        )
    }

    pub fn get_column(&'_ self, column: usize) -> ViewBigVector<'_> {
        ViewBigVector::new(self.row_count, &self.numbers, column, self.column_count)
    }

    pub fn get_column_mut(&'_ mut self, column: usize) -> MutViewBigVector<'_> {
        MutViewBigVector::new(self.row_count, &mut self.numbers, column, self.column_count)
    }

    pub fn set_row(&mut self, row: usize, value: OwnedBigVector) {
        let mut i = 0;
        for x in value.numbers {
            self.set(row, i, x);
            i += 1;
        }
    }

    pub fn set_column(&mut self, column: usize, value: OwnedBigVector) {
        let mut i = 0;
        for x in value.numbers {
            self.set(i, column, x);
            i += 1;
        }
    }

    pub fn shink(
        &mut self,
        start_row: usize,
        start_column: usize,
        row_count: usize,
        column_count: usize,
    ) {
        if start_row == 0
            && start_column == 0
            && self.row_count == row_count
            && self.column_count == column_count
        {
            return;
        }

        let old_numbers = std::mem::replace(
            &mut self.numbers,
            Vec::with_capacity(row_count * column_count),
        );

        let start = start_row * self.column_count + start_column;
        let skip = self.column_count - column_count;

        let mut iter = old_numbers.into_iter().skip(start);

        for _ in 0..row_count {
            for _ in 0..column_count {
                self.numbers.push(iter.next().unwrap());
            }

            for _ in 0..skip {
                iter.next();
            }
        }

        self.row_count = row_count;
        self.column_count = column_count;
    }

    pub fn transpose(&self) -> BigMatrix {
        let mut dest: BigMatrix = BigMatrix::new(self.row_count, self.column_count);

        for i in 0..self.column_count {
            dest.set_row(i, self.get_column(i).to_owned());
        }

        dest
    }

    pub fn swap_rows(&mut self, row1: usize, row2: usize) {
        let n = &mut self.numbers;

        for i in 0..self.column_count {
            n.swap(i + self.column_count * row1, i + self.column_count * row2);
        }
    }

    pub fn swap_elements(&mut self, row1: usize, col1: usize, row2: usize, col2: usize) {
        self.numbers.swap(
            col1 + self.column_count * row1,
            col2 + self.column_count * row2,
        );
    }
}

macro_rules! impl_for_owned_and_refs {
    ($macro: ident, $a: ty) => {
        $macro!($a, BigMatrix);
        $macro!(&$a, BigMatrix);
        $macro!($a, &BigMatrix);
        $macro!(&$a, &BigMatrix);
    };
}

macro_rules! impl_operation {
    ($a: ty, $b: ty, $t: ident, $m: ident, $f: tt) => {
        impl $t<$a> for $b {

            type Output = BigMatrix;

            fn $m(self, other: $a) -> BigMatrix {
                let mut vec: Vec<Rational> = Vec::with_capacity(self.row_count * self.column_count);
                for row in 0..self.row_count {
                    for col in 0..self.column_count {
                        vec.push(&*self.get(row, col) $f other.get(row, col))
                    }
                }

                BigMatrix {
                    numbers: vec,
                    row_count: self.row_count,
                    column_count: self.column_count
                }
            }

        }
    };
}

macro_rules! add {
    ($a: ty, $b: ty) => {
        impl_operation!($a, $b, Add, add, +);
    };
}
impl_for_owned_and_refs!(add, BigMatrix);

macro_rules! sub {
    ($a: ty, $b: ty) => {
        impl_operation!($a, $b, Sub, sub, -);
    };
}
impl_for_owned_and_refs!(sub, BigMatrix);

macro_rules! mul_with_vec {
    ($a: ty, $b: ty) => {
        impl Mul<$a> for $b {
            type Output = OwnedBigVector;

            fn mul(self, v: $a) -> OwnedBigVector {
                OwnedBigVector::new_from(self.row_count, |i| v.dot(&self.get_row(i)))
            }
        }
    };
}

impl_for_owned_and_refs!(mul_with_vec, OwnedBigVector);

macro_rules! mul_with_matrix {
    ($a: ty, $b: ty) => {
        impl Mul<$a> for $b {
            type Output = BigMatrix;

            fn mul(self, m: $a) -> BigMatrix {
                BigMatrix::new_provided(self.row_count, self.column_count, |row, col| {
                    self.get_row(row).dot(&m.get_column(col))
                })
            }
        }
    };
}

impl_for_owned_and_refs!(mul_with_matrix, BigMatrix);

impl Clone for BigMatrix {
    fn clone(&self) -> BigMatrix {
        BigMatrix {
            numbers: self.numbers.clone(),
            row_count: self.row_count,
            column_count: self.column_count,
        }
    }
}

impl PartialEq for BigMatrix {
    fn eq(&self, other: &Self) -> bool {
        if self.row_count != other.row_count || self.column_count != other.column_count {
            return false;
        }

        for row in 0..self.row_count {
            for col in 0..self.column_count {
                if self.get(row, col) != other.get(row, col) {
                    return false;
                }
            }
        }

        true
    }
}
