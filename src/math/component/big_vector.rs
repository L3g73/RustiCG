use malachite::base::num::basic::traits::{One, Zero};
use malachite::rational::Rational;
use std::ops::AddAssign;
use std::ops::DivAssign;
use std::ops::MulAssign;
use std::ops::SubAssign;
use std::ops::{Add, Mul, Sub};

pub trait Vector {
    fn get(&self, index: usize) -> &Rational;
    fn dimension(&self) -> usize;
}

impl<T: Vector> Vector for &T {
    fn get(&self, index: usize) -> &Rational {
        T::get(self, index)
    }

    fn dimension(&self) -> usize {
        T::dimension(self)
    }
}

pub trait Dot<T> {
    fn dot(&self, other: &T) -> Rational;
}

impl<T: Vector, S: Vector> Dot<T> for S {
    fn dot(&self, other: &T) -> Rational {
        let mut result = Rational::ZERO;

        for i in 0..self.dimension() {
            result += self.get(i) * other.get(i);
        }

        result
    }
}

pub trait MagnitudeSq {
    fn magnitude_sq(&self) -> Rational;
}

impl<T: Dot<T>> MagnitudeSq for T {
    fn magnitude_sq(&self) -> Rational {
        self.dot(self)
    }
}

pub trait IsZero {
    fn is_zero(&self) -> bool;
}

impl<T: Vector> IsZero for T {
    fn is_zero(&self) -> bool {
        for i in 0..self.dimension() {
            if *self.get(i) != 0 {
                return false;
            }
        }

        true
    }
}

macro_rules! impl_vec_partial_eq {
    (
        $name:ident
        $(<$( $lt:lifetime ),+>)?
    ) => {
        impl<$($( $lt ),+ ,)? T: Vector> PartialEq<T> for $name $(< $( $lt ),+ >)? {
            fn eq(&self, other: &T) -> bool {
                if self.dimension() != other.dimension() {
                    return false;
                }

                for i in 0..self.dimension() {
                    if self.get(i) != other.get(i) {
                        return false;
                    }
                }

                true
            }
        }
    };
}

macro_rules! impl_vec_traits {
    (
        $name:ident
        $(<$( $lt:lifetime ),+>)?
    ) => {

        macro_rules! combine_vec {
            ($trait_name: ident, $method: ident) => {
                impl<$($( $lt ),+ ,)? T: Vector> $trait_name<T> for $name $(< $( $lt ),+ >)? {
                    type Output = OwnedBigVector;
                    fn $method(self, other: T) -> OwnedBigVector {
                        OwnedBigVector::new_from(self.dimension, |i| { self.get(i).$method(other.get(i)) })
                    }
                }
                impl<$($( $lt ),+ ,)? T: Vector> $trait_name<T> for &$name $(< $( $lt ),+ >)? {
                    type Output = OwnedBigVector;
                    fn $method(self, other: T) -> OwnedBigVector {
                        OwnedBigVector::new_from(self.dimension, |i| { self.get(i).$method(other.get(i)) })
                    }
                }
            };
        }

        combine_vec!(Add, add);
        combine_vec!(Sub, sub);

        impl$(< $( $lt ),+ >)? Mul<&Rational> for $name $(< $( $lt ),+ >)? {
            type Output = OwnedBigVector;

            fn mul(self, other: &Rational) -> OwnedBigVector {
                OwnedBigVector::new_from(self.dimension, |i| {
                    self.get(i) * other
                })
            }
        }
        impl$(< $( $lt ),+ >)? Mul<&Rational> for &$name $(< $( $lt ),+ >)? {
            type Output = OwnedBigVector;

            fn mul(self, other: &Rational) -> OwnedBigVector {
                OwnedBigVector::new_from(self.dimension, |i| {
                    self.get(i) * other
                })
            }
        }
    };
}

macro_rules! impl_mut_vec_traits {
    (
        $name:ident
        $(<$( $lt:lifetime ),+>)?
    ) => {
        macro_rules! assign_vec {
            ($ty: ident, $method: ident) => {
                impl<$($( $lt ),+ ,)? T: Vector> $ty<T> for $name $(< $( $lt ),+ >)? {
                    fn $method(&mut self, other: T) {
                        for i in 0..self.dimension {
                            self.get_mut(i).$method(other.get(i));
                        }
                    }
                }
            };
        }
        assign_vec!(AddAssign, add_assign);
        assign_vec!(SubAssign, sub_assign);

        macro_rules! assign_scalar {
            ($ty: ident, $method: ident) => {
                impl $(< $( $lt ),+ >)? $ty<&Rational> for $name $(< $( $lt ),+ >)? {

                    fn $method(&mut self, other: &Rational) {
                        for i in 0..self.dimension {
                            self.get_mut(i).$method(other);
                        }
                    }

                }
            };
        }
        assign_scalar!(MulAssign, mul_assign);
        assign_scalar!(DivAssign, div_assign);
    };
}

pub struct OwnedBigVector {
    pub dimension: usize,
    pub(crate) numbers: Vec<Rational>,
}

impl Vector for OwnedBigVector {
    fn get(&self, index: usize) -> &Rational {
        &self.numbers[index]
    }

    fn dimension(&self) -> usize {
        self.dimension
    }
}

impl OwnedBigVector {
    pub fn new_from<F: Fn(usize) -> Rational>(dimension: usize, f: F) -> OwnedBigVector {
        let mut vec: Vec<Rational> = Vec::with_capacity(dimension);

        for i in 0..dimension {
            vec.push(f(i));
        }

        Self::new_from_numbers(vec)
    }

    pub fn new(dimension: usize) -> OwnedBigVector {
        Self::new_from(dimension, |_i| Rational::ZERO)
    }

    pub fn new_from_numbers(numbers: Vec<Rational>) -> OwnedBigVector {
        OwnedBigVector {
            dimension: numbers.len(),
            numbers,
        }
    }

    pub fn as_view(&'_ self) -> ViewBigVector<'_> {
        ViewBigVector::new(self.dimension, &self.numbers, 0, 1)
    }

    pub fn basis(size: usize, i: usize) -> OwnedBigVector {
        Self::basis_scaled(size, i, Rational::ONE)
    }

    pub fn basis_scaled(size: usize, i: usize, scale: Rational) -> OwnedBigVector {
        let mut vector = Self::new(size);
        vector.set(i, scale);
        vector
    }

    fn get_mut(&mut self, index: usize) -> &mut Rational {
        &mut self.numbers[index]
    }

    pub fn set(&mut self, index: usize, value: Rational) {
        self.numbers[index] = value
    }
}

impl_vec_partial_eq!(OwnedBigVector);
impl_mut_vec_traits!(OwnedBigVector);

macro_rules! combine_vec {
    ($trait_name: ident, $method: ident, $method_assign: ident) => {
        impl<T: Vector> $trait_name<T> for OwnedBigVector {
            type Output = OwnedBigVector;
            fn $method(mut self, other: T) -> OwnedBigVector {
                self.$method_assign(other);
                self
            }
        }
        impl<T: Vector> $trait_name<T> for &OwnedBigVector {
            type Output = OwnedBigVector;
            fn $method(self, other: T) -> OwnedBigVector {
                OwnedBigVector::new_from(self.dimension, |i| self.get(i).$method(other.get(i)))
            }
        }
    };
}

combine_vec!(Add, add, add_assign);
combine_vec!(Sub, sub, sub_assign);

impl Mul<&Rational> for OwnedBigVector {
    type Output = OwnedBigVector;

    fn mul(mut self, other: &Rational) -> OwnedBigVector {
        self.mul_assign(other);
        self
    }
}
impl Mul<&Rational> for &OwnedBigVector {
    type Output = OwnedBigVector;

    fn mul(self, other: &Rational) -> OwnedBigVector {
        OwnedBigVector::new_from(self.dimension, |i| self.get(i) * other)
    }
}

impl Clone for OwnedBigVector {
    fn clone(&self) -> Self {
        OwnedBigVector {
            dimension: self.dimension,
            numbers: self.numbers.clone(),
        }
    }
}

pub struct ViewBigVector<'r> {
    pub dimension: usize,
    pub numbers: &'r Vec<Rational>,
    start_pos: usize,
    step: usize,
}

impl Vector for ViewBigVector<'_> {
    fn get(&self, index: usize) -> &Rational {
        &self.numbers[self.step * index + self.start_pos]
    }

    fn dimension(&self) -> usize {
        self.dimension
    }
}

impl<'r> ViewBigVector<'r> {
    pub fn new(
        dimension: usize,
        numbers: &'r Vec<Rational>,
        start_pos: usize,
        step: usize,
    ) -> ViewBigVector<'r> {
        ViewBigVector {
            dimension,
            numbers,
            start_pos,
            step,
        }
    }

    pub fn to_owned(&self) -> OwnedBigVector {
        let mut vec: Vec<Rational> = Vec::with_capacity(self.dimension);
        for i in 0..self.dimension {
            vec.push(self.get(i).clone())
        }

        OwnedBigVector::new_from_numbers(vec)
    }
}

impl_vec_partial_eq!(ViewBigVector<'r>);
impl_vec_traits!(ViewBigVector<'r>);

pub struct MutViewBigVector<'r> {
    pub dimension: usize,
    pub numbers: &'r mut Vec<Rational>,
    start_pos: usize,
    step: usize,
}

impl Vector for MutViewBigVector<'_> {
    fn get(&self, index: usize) -> &Rational {
        &self.numbers[self.step * index + self.start_pos]
    }

    fn dimension(&self) -> usize {
        self.dimension
    }
}

impl<'r> MutViewBigVector<'r> {
    pub fn new(
        dimension: usize,
        numbers: &'r mut Vec<Rational>,
        start_pos: usize,
        step: usize,
    ) -> MutViewBigVector<'r> {
        MutViewBigVector {
            dimension,
            numbers,
            start_pos,
            step,
        }
    }

    fn get_mut(&mut self, index: usize) -> &mut Rational {
        &mut self.numbers[self.step * index + self.start_pos]
    }

    pub fn set(&mut self, index: usize, value: Rational) {
        self.numbers[self.step * index + self.start_pos] = value
    }
}

impl_vec_partial_eq!(MutViewBigVector<'r>);
impl_vec_traits!(MutViewBigVector<'r>);
impl_mut_vec_traits!(MutViewBigVector<'r>);
