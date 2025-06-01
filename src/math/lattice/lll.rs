use crate::math::component::big_matrix::BigMatrix;
use crate::math::component::big_vector::OwnedBigVector;
use crate::math::component::big_vector::{Dot, IsZero, MagnitudeSq, Vector};
use malachite::Integer;
use malachite::base::num::arithmetic::traits::Floor;
use malachite::base::num::basic::traits::{One, OneHalf, Zero};
use malachite::rational::Rational;
use std::cmp::max;
use std::ops::{MulAssign, SubAssign};

pub struct Params {
    pub max_stage: i32,
    pub delta: Rational,
}

pub struct LLL {
    nb_rows: usize,
    nb_cols: usize,
    coordinates: BigMatrix,
    base_gso: BigMatrix,
    mu: BigMatrix,
    norms: OwnedBigVector,
    basis: BigMatrix,
    params: Params,
}

pub fn reduce(lattice: BigMatrix) -> BigMatrix {
    LLL {
        nb_rows: lattice.row_count,
        nb_cols: lattice.column_count,
        base_gso: BigMatrix::new(lattice.row_count, lattice.column_count),
        mu: BigMatrix::new(lattice.row_count, lattice.row_count),
        norms: OwnedBigVector::new(lattice.row_count),
        coordinates: BigMatrix::identity(lattice.row_count),
        basis: lattice,
        params: Params {
            max_stage: -1,
            delta: Rational::from_signeds(99, 100),
        },
    }
    .reduce_lll()
}

impl LLL {
    fn test_condition(&self, k: usize) -> bool {
        let mu_temp = self.mu.get(k, k - 1);
        let factor = &self.params.delta - mu_temp * mu_temp;
        self.norms.get(k) < &(self.norms.get(k - 1) * factor)
    }

    fn update_gso(&mut self, k: usize) {
        let mut new_row = self.basis.get_row(k).to_owned();

        for j in 0..k {
            if self.norms.get(j) != &Rational::ZERO {
                self.mu.set(
                    k,
                    j,
                    self.basis.get_row(k).dot(&self.base_gso.get_row(j)) / self.norms.get(j),
                )
            } else {
                self.mu.set(k, j, Rational::ZERO);
            }

            new_row -= self.base_gso.get_row(j) * self.mu.get(k, j);
        }

        self.norms.set(k, new_row.magnitude_sq());
        self.base_gso.set_row(k, new_row);
    }

    fn red(&mut self, i: usize, j: usize) {
        let r = &(self.mu.get(i, j) + Rational::ONE_HALF).floor();
        if r == &Integer::ZERO {
            return;
        }

        let r = &Rational::from(r);

        let scaled_basis_row = self.basis.get_row(j) * r;
        self.basis.get_row_mut(i).sub_assign(scaled_basis_row);
        let scaled_coordinates_row = self.coordinates.get_row(j) * r;
        self.coordinates
            .get_row_mut(i)
            .sub_assign(scaled_coordinates_row);

        self.mu.set(i, j, self.mu.get(i, j) - r);

        for col in 0..j {
            self.mu
                .set(i, col, self.mu.get(i, col) - self.mu.get(j, col) * r);
        }
    }

    fn swapg(&mut self, k: usize, kmax: usize) {
        self.basis.swap_rows(k, k - 1);
        self.coordinates.swap_rows(k, k - 1);

        if k > 1 {
            for j in 0..k - 1 {
                self.mu.swap_elements(k, j, k - 1, j);
            }
        }

        let tmu = &self.mu.get(k, k - 1).clone();
        let t_b = self.norms.get(k) + tmu * tmu * self.norms.get(k - 1);

        if t_b == Rational::ZERO {
            self.norms.set(k, self.norms.get(k - 1).clone());
            self.norms.set(k - 1, Rational::ZERO);
            self.base_gso.swap_rows(k, k - 1);

            for i in k + 1..=kmax {
                self.mu.set(i, k, self.mu.get(i, k - 1).clone());
                self.mu.set(i, k - 1, Rational::ZERO);
            }
        } else if self.norms.get(k) == &Rational::ZERO && tmu != &Rational::ZERO {
            self.norms.set(k - 1, t_b);
            self.base_gso.get_row_mut(k - 1).mul_assign(tmu);
            self.mu.set(k, k - 1, Rational::ONE / tmu);

            for i in k + 1..=kmax {
                self.mu.set(i, k, self.mu.get(i, k - 1) / tmu);
            }
        } else {
            let t = &(self.norms.get(k - 1) / &t_b);
            self.mu.set(k, k - 1, tmu * t);
            let b = self.base_gso.get_row(k - 1).to_owned();
            self.base_gso
                .set_row(k - 1, self.base_gso.get_row(k) + &b * tmu);
            self.base_gso.set_row(
                k,
                b * &(self.norms.get(k) / &t_b) - self.base_gso.get_row(k) * self.mu.get(k, k - 1),
            );

            self.norms.set(k, self.norms.get(k) * t);
            self.norms.set(k - 1, t_b);

            for i in k + 1..=kmax {
                let new_t = &self.mu.get(i, k).clone();
                self.mu.set(i, k, self.mu.get(i, k - 1) - (tmu * new_t));
                self.mu
                    .set(i, k - 1, self.mu.get(k, k - 1) * self.mu.get(i, k) + new_t);
            }
        }
    }

    fn remove_zeroes(&mut self) {
        let mut p: usize = 0;

        for i in 0..self.nb_rows {
            if self.basis.get_row(i).is_zero() {
                p += 1;
            }
        }

        self.basis.shink(p, 0, self.nb_rows - p, self.nb_cols);
        self.base_gso.shink(p, 0, self.nb_rows - p, self.nb_cols);
        self.mu.shink(p, p, self.nb_rows - p, self.nb_rows - p);

        self.norms = OwnedBigVector::new_from(self.nb_rows - p, |i| self.norms.get(i + p).clone());
    }

    fn reduce_lll(mut self) -> BigMatrix {
        let mut k = 1;
        let mut kmax = 0;
        let n = if self.params.max_stage == -1 {
            self.nb_rows
        } else {
            self.params.max_stage as usize
        };

        self.base_gso.set_row(0, self.basis.get_row(0).to_owned());
        let mut update_gso = true;
        self.norms.set(0, self.basis.get_row(0).magnitude_sq());

        while k < n {
            if k > kmax && update_gso {
                kmax = k;
                self.update_gso(k);
            }

            self.red(k, k - 1);
            if self.test_condition(k) {
                self.swapg(k, kmax);
                k = max(1, k - 1);
                update_gso = false;
            } else {
                let mut l = k as i32 - 2;
                while l >= 0 {
                    self.red(k, l as usize);
                    l -= 1;
                }
                k += 1;
                update_gso = true;
            }
        }

        self.remove_zeroes();
        self.basis
    }
}
