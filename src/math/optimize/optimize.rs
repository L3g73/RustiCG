use crate::math::component::big_matrix::BigMatrix;
use crate::math::component::big_vector::{OwnedBigVector, Vector, ViewBigVector};
use crate::math::optimize::gauss_jordan::{reduce_matrix, reduce_with_action};
use crate::util::util::{copy_of_usizes, filled_vec};
use malachite::base::num::arithmetic::traits::{Reciprocal, Sign};
use malachite::base::num::basic::traits::{NegativeOne, One, Zero};
use malachite::rational::Rational;
use std::cmp::Ordering::{Greater, Less};
use std::ops::{AddAssign, DivAssign, Mul, MulAssign, Neg, SubAssign};
use std::sync::Arc;

pub struct Optimize {
    transform: Arc<BigMatrix>,
    table: BigMatrix,
    basics: Vec<usize>,
    nonbasics: Vec<usize>,

    rows: usize,
    rows_m1: usize,
    cols: usize,
    cols_m1: usize,
}

impl Optimize {
    const NO_VALUE: usize = usize::MAX;

    pub fn new(
        table: BigMatrix,
        basics: Vec<usize>,
        nonbasics: Vec<usize>,
        transform: Arc<BigMatrix>,
    ) -> Optimize {
        let rows = table.row_count;
        let cols = table.column_count;
        Optimize {
            transform,
            table,
            basics,
            nonbasics,
            rows,
            rows_m1: rows - 1,
            cols,
            cols_m1: cols - 1,
        }
    }

    fn transform_for_table(&self, lhs: &ViewBigVector, rhs: Rational) -> OwnedBigVector {
        let mut transformed = OwnedBigVector::new(self.transform.column_count);
        let mut eliminated = OwnedBigVector::new(self.cols);

        transformed.set(self.transform.column_count - 1, rhs);

        for row in 0..self.transform.row_count {
            let x = lhs.get(row);
            transformed -= self.transform.get_row(row) * x;
        }

        for col in 0..self.cols_m1 {
            eliminated.set(col, transformed.get(self.nonbasics[col]).clone());
        }

        eliminated.set(
            self.cols_m1,
            transformed.get(self.transform.column_count - 1).clone(),
        );

        for row in 0..self.rows_m1 {
            let x = transformed.get(self.basics[row]);
            eliminated -= self.table.get_row(row) * x;
        }

        eliminated
    }

    pub fn maximize(&mut self, gradient: &ViewBigVector) -> Rational {
        let result = self.minimize(&(gradient * &Rational::NEGATIVE_ONE).as_view());
        result.neg()
    }

    pub fn minimize(&mut self, gradient: &ViewBigVector) -> Rational {
        let mut temp = OwnedBigVector::new(self.cols);
        temp -= &self.transform_for_table(gradient, Rational::ZERO);
        self.table.set_row(self.rows_m1, temp);

        self.solve();

        self.table.get(self.rows_m1, self.cols_m1).clone()
    }

    fn solve(&mut self) {
        loop {
            if !self.step() {
                break;
            }
        }
    }

    fn step(&mut self) -> bool {
        let mut bland = false;

        let mut entering: i32 = -1;
        let mut exiting: i32 = -1;

        let mut candidate = Rational::ZERO;

        for row in 0..self.rows_m1 {
            if *self.table.get(row, self.cols_m1) == 0 {
                bland = true;
                break;
            }
        }

        for col in 0..self.cols_m1 {
            let x = self.table.get(self.rows_m1, col);
            if x.sign() != Greater || (entering != -1 && *x <= candidate) {
                continue;
            }

            entering = col as i32;
            candidate = x.clone();

            if bland {
                break;
            }
        }

        if entering == -1 {
            return false;
        }

        for row in 0..self.rows_m1 {
            let x = self.table.get(row, entering as usize);
            if x.sign() != Greater {
                continue;
            }

            let y = self.table.get(row, self.cols_m1) / x;

            if exiting != -1 && y >= candidate {
                continue;
            }

            exiting = row as i32;
            candidate = y;
        }

        self.pivot(entering as usize, exiting as usize);
        true
    }

    fn pivot(&mut self, entering: usize, exiting: usize) {
        let pivot = self.table.get(exiting, entering).clone();

        for col in 0..self.cols {
            if col == entering {
                continue;
            }

            self.table
                .set(exiting, col, self.table.get(exiting, col) / &pivot);
        }

        for row in 0..self.rows {
            if row == exiting {
                continue;
            }

            let x = self.table.get(row, entering).clone();

            for col in 0..self.cols {
                if col == entering {
                    continue;
                }

                let y = self.table.get(exiting, col);
                self.table
                    .set(row, col, self.table.get(row, col) - (&x * y))
            }

            self.table.set(row, entering, -&(x / &pivot));
        }

        self.table.set(exiting, entering, pivot.reciprocal());

        let new_basic = self.nonbasics[entering];
        self.nonbasics[entering] = self.basics[exiting];
        self.basics[exiting] = new_basic;
    }

    pub fn with_strict_bound(&self, lhs: &ViewBigVector, rhs: Rational) -> Optimize {
        let mut new_table = BigMatrix::new(self.rows + 1, self.cols);

        for row in 0..self.rows_m1 {
            new_table.set_row(row, self.table.get_row(row).to_owned());
        }

        new_table.set_row(self.rows_m1, self.transform_for_table(lhs, rhs));

        if new_table.get(self.rows_m1, self.cols_m1).sign() == Less {
            new_table
                .get_row_mut(self.rows_m1)
                .mul_assign(&Rational::NEGATIVE_ONE);
        }

        let mut new_basics = copy_of_usizes(&self.basics, self.rows);
        let new_nonbasics = copy_of_usizes(&self.nonbasics, self.cols_m1);

        new_basics[self.rows_m1] = self.rows_m1 + self.cols_m1;

        Self::from(
            new_table,
            new_basics,
            new_nonbasics,
            1,
            self.transform.clone(),
        )
    }

    pub fn from(
        mut table: BigMatrix,
        basics: Vec<usize>,
        nonbasics: Vec<usize>,
        artificials: usize,
        transform: Arc<BigMatrix>,
    ) -> Optimize {
        let rows = table.row_count;
        let rows_m1 = rows - 1;
        let cols = table.column_count;
        let cols_m1 = cols - 1;

        let real_variables = rows_m1 + cols_m1 - artificials;

        for basic_row in 0..rows_m1 {
            if basics[basic_row] < real_variables {
                continue;
            }

            let row = table.get_row(basic_row).to_owned();
            table.get_row_mut(rows_m1).add_assign(row);
        }

        let mut optimize = Optimize::new(table, basics, nonbasics, Arc::new(BigMatrix::new(0, 0)));
        optimize.solve();

        for row in 0..rows_m1 {
            if optimize.basics[row] >= real_variables {
                for col in 0..cols_m1 {
                    if optimize.nonbasics[col] >= real_variables
                        || *optimize.table.get(row, col) == 0
                    {
                        continue;
                    }

                    optimize.pivot(col, row);
                    break;
                }
            }
        }

        let final_cols = cols - artificials;
        let mut final_table = BigMatrix::new(rows, final_cols);

        let mut c0: usize = 0;
        let mut c1: usize = 0;

        while c0 < final_cols - 1 {
            loop {
                if optimize.nonbasics[c1] >= real_variables {
                    c1 += 1;
                    continue;
                }

                for row in 0..rows_m1 {
                    final_table.set(row, c0, optimize.table.get(row, c1).clone());
                }

                optimize.nonbasics[c0] = optimize.nonbasics[c1];
                break;
            }

            c0 += 1;
            c1 += 1;
        }

        for row in 0..rows_m1 {
            final_table.set(
                row,
                final_cols - 1,
                optimize.table.get(row, cols_m1).clone(),
            );
        }

        Optimize::new(final_table, optimize.basics, optimize.nonbasics, transform)
    }

    pub fn from_build(mut inner_table: BigMatrix, transform: Arc<BigMatrix>) -> Optimize {
        let constraints = inner_table.row_count;
        let variables = inner_table.column_count - 1;

        let mut basics: Vec<usize> = filled_vec(constraints, |_i| Self::NO_VALUE);

        let mut nonbascis: Vec<usize> = Vec::new();

        for row in 0..constraints {
            if inner_table.get(row, variables).sign() != Less {
                continue;
            }

            inner_table
                .get_row_mut(row)
                .mul_assign(&Rational::NEGATIVE_ONE);
        }

        for col in 0..variables {
            let mut count = 0;
            let mut index = Self::NO_VALUE;

            for row in 0..inner_table.row_count {
                if *inner_table.get(row, col) == 0 {
                    continue;
                }

                count += 1;
                index = row;
            }

            if count == 1
                && basics[index] == Self::NO_VALUE
                && inner_table.get(index, col).sign() == Greater
            {
                let divisor = inner_table.get(index, col).clone();
                inner_table.get_row_mut(index).div_assign(&divisor);
                basics[index] = col;
            } else {
                nonbascis.push(col);
            }
        }

        let mut artificials = 0;

        for row in 0..constraints {
            if basics[row] != Self::NO_VALUE {
                continue;
            }

            basics[row] = variables + artificials;
            artificials += 1;
        }

        let nonbasic_count = variables - constraints + artificials;
        let mut table = BigMatrix::new(constraints + 1, nonbasic_count + 1);

        for row in 0..constraints {
            for basic_row in 0..constraints {
                if row == basic_row || basics[basic_row] >= variables {
                    continue;
                }

                let scalar = inner_table.get(row, basics[basic_row]);
                let basic_vector = inner_table.get_row(basic_row).mul(scalar);
                let mut row_vector = inner_table.get_row_mut(row);

                row_vector.sub_assign(basic_vector.mul(row_vector.get(basics[basic_row])))
            }

            for col in 0..nonbasic_count {
                table.set(row, col, inner_table.get(row, nonbascis[col]).clone());
            }

            table.set(row, nonbasic_count, inner_table.get(row, variables).clone());
        }

        Self::from(table, basics, nonbascis, artificials, transform)
    }
}

impl Clone for Optimize {
    fn clone(&self) -> Self {
        Optimize::new(
            self.table.clone(),
            copy_of_usizes(&self.basics, self.rows_m1),
            copy_of_usizes(&self.nonbasics, self.cols_m1),
            self.transform.clone(),
        )
    }
}

pub struct OptimizeBuilder {
    size: usize,
    slacks: Vec<i32>,
    lefts: Vec<OwnedBigVector>,
    rights: Vec<Rational>,
}

impl OptimizeBuilder {
    pub fn of_size(size: usize) -> OptimizeBuilder {
        OptimizeBuilder {
            size,
            slacks: Vec::new(),
            lefts: Vec::new(),
            rights: Vec::new(),
        }
    }

    fn add(&mut self, slack: i32, lhs: OwnedBigVector, rhs: Rational) {
        self.slacks.push(slack);
        self.lefts.push(lhs);
        self.rights.push(rhs);
    }

    pub fn with_lower_bound(&mut self, lhs: usize, rhs: Rational) {
        self.add(-1, OwnedBigVector::basis(self.size, lhs), rhs);
    }

    pub fn with_upper_bound(&mut self, lhs: usize, rhs: Rational) {
        self.add(1, OwnedBigVector::basis(self.size, lhs), rhs);
    }

    pub fn build(self) -> Optimize {
        let column_count = (self.size + self.slacks.len()) * 2 + self.size;
        let mut constraint = 0;
        let mut slack = self.size;

        let mut table = BigMatrix::new(self.slacks.len() + self.size, column_count + 1);

        while constraint < self.slacks.len() {
            for col in 0..self.size {
                table.set(constraint, col, self.lefts[constraint].get(col).clone())
            }

            table.set(constraint, column_count, self.rights[constraint].clone());

            if self.slacks[constraint] != 0 {
                table.set(constraint, slack, Rational::from(self.slacks[constraint]));
                slack += 1;
            }

            constraint += 1;
        }

        let mut privot_rows =
            reduce_with_action(&mut table, |col: usize, _rows: &Vec<i32>| col < self.size);

        for col in 0..self.size {
            if privot_rows[col] != -1 {
                continue;
            }

            let mut row = table.get_row_mut(constraint);
            row.set(col, Rational::ONE);
            row.set(slack, Rational::ONE);
            row.set(slack + 1, Rational::NEGATIVE_ONE);

            constraint += 1;
            slack += 1;
        }

        privot_rows = reduce_matrix(&mut table);

        constraint = (1 + privot_rows.iter().max().or_else(|| Some(&-1)).unwrap()) as usize;

        let mut transform = BigMatrix::new(self.size, slack - self.size + 1);
        let mut inner_table = BigMatrix::new(constraint - self.size, slack - self.size + 1);

        for row in 0..self.size {
            for col in 0..slack - self.size {
                transform.set(row, col, table.get(row, self.size + col).clone());
            }

            transform.set(row, slack - self.size, table.get(row, column_count).clone());
        }

        for row in 0..self.size {
            for col in 0..slack - self.size {
                inner_table.set(
                    row,
                    col,
                    table.get(self.size + row, self.size + col).clone(),
                );
            }

            inner_table.set(
                row,
                slack - self.size,
                table.get(self.size + row, column_count).clone(),
            );
        }

        Optimize::from_build(inner_table, Arc::new(transform))
    }
}
