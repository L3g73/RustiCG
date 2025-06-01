use malachite::base::num::arithmetic::traits::ExtendedGcd;
use malachite::base::num::basic::traits::Zero;
use malachite::{Integer, Natural};

pub fn filled_vec<T, F>(capacity: usize, value_func: F) -> Vec<T>
where
    F: Fn(usize) -> T,
{
    let mut vec = Vec::with_capacity(capacity);
    for i in 0..capacity {
        vec.push(value_func(i));
    }

    vec
}

pub fn copy_of_usizes(vec: &Vec<usize>, new_size: usize) -> Vec<usize> {
    let mut result = Vec::with_capacity(new_size);
    let old_size = vec.len();

    result.extend_from_slice(&vec[0..new_size.min(old_size)]);

    if new_size > old_size {
        result.resize_with(new_size, Default::default);
    }

    result
}

pub trait Ushr {
    fn ushr(self, n: i32) -> Self;
}

impl Ushr for i64 {
    fn ushr(self, n: i32) -> Self {
        if self >= 0 {
            return self >> n;
        }

        let mut u = u64::from_ne_bytes(self.to_ne_bytes());
        u >>= n;
        i64::from_ne_bytes(u.to_ne_bytes())
    }
}

pub trait LCM {
    fn lcm(&self, other: &Self) -> Self;
}

impl LCM for Integer {
    fn lcm(&self, other: &Self) -> Self {
        (self * other) / Integer::from(self.extended_gcd(other).0)
    }
}

pub trait ToU64 {
    fn to_u64(&self) -> u64;
}

impl ToU64 for Natural {
    fn to_u64(&self) -> u64 {
        if self == &Natural::ZERO {
            return 0_u64;
        }

        self.to_limbs_asc()[0]
    }
}
