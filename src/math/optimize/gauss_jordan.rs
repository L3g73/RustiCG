use crate::math::component::big_matrix::BigMatrix;
use std::ops::{DivAssign, SubAssign};

#[inline(always)]
fn for_all<T: Fn(&mut BigMatrix)>(matrix: &mut BigMatrix, other: &mut Vec<BigMatrix>, action: T) {
    action(matrix);
    for x in other {
        action(x);
    }
}

pub fn reduce_with_action<T: Fn(usize, &Vec<i32>) -> bool>(
    matrix: &mut BigMatrix,
    action: T,
) -> Vec<i32> {
    reduce(matrix, &mut Vec::new(), action)
}

pub fn reduce_matrix(matrix: &mut BigMatrix) -> Vec<i32> {
    reduce(matrix, &mut Vec::new(), |_col: usize, _rows: &Vec<i32>| {
        true
    })
}

pub fn reduce<T: Fn(usize, &Vec<i32>) -> bool>(
    matrix: &mut BigMatrix,
    others: &mut Vec<BigMatrix>,
    action: T,
) -> Vec<i32> {
    let mut pivot_rows: Vec<i32> = Vec::with_capacity(matrix.column_count);
    for _i in 0..matrix.column_count {
        pivot_rows.push(-1);
    }

    let mut row: usize = 0;
    let mut pivot_column: usize = 0;

    while row < matrix.row_count && pivot_column < matrix.column_count {
        let mut pivot_row = row;
        while pivot_row < matrix.row_count {
            if *matrix.get(pivot_row, pivot_column) != 0 {
                break;
            }

            pivot_row += 1;
        }

        if pivot_row < matrix.row_count {
            let pivot = matrix.get(pivot_row, pivot_column).clone();

            for_all(matrix, others, |matrix| {
                matrix.get_row_mut(pivot_row).div_assign(&pivot)
            });

            for i in 0..matrix.row_count {
                if i == pivot_row {
                    continue;
                }

                let scale = matrix.get(i, pivot_column).clone();
                for_all(matrix, others, |matrix| {
                    let scaled_row = matrix.get_row(pivot_row) * &scale;
                    matrix.get_row_mut(i).sub_assign(scaled_row);
                });
            }

            for_all(matrix, others, |matrix| matrix.swap_rows(row, pivot_row));
            pivot_rows[pivot_column] = row as i32;
            row += 1;
        }

        loop {
            pivot_column += 1;
            if pivot_column >= matrix.column_count || action(pivot_column, &pivot_rows) {
                break;
            }
        }
    }

    pivot_rows
}
