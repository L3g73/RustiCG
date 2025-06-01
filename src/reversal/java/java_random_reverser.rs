use crate::reversal::random_reverser::RandomReverser;
use crate::util::util::Ushr;

pub fn add_bound_next_int_call(
    reverser: &mut RandomReverser,
    bound: i32,
    min: Vec<i32>,
    max: Vec<i32>,
) {
    if (bound & -bound) == bound {
        // bound is power of 2
        let log = bound.trailing_zeros() as i64;
        let offset = 1_i64 << (48 - log);

        let min: Vec<i64> = min.into_iter().map(|min| (min as i64) * offset).collect();
        let max: Vec<i64> = max.into_iter().map(|max| (max as i64) * offset + offset - 1).collect();
        reverser.add_measured_seed(min, max);
    } else {
        let offset = 1_i64 << 17;

        let min: Vec<i64> = min.into_iter().map(|min| (min as i64) * offset).collect();
        let max: Vec<i64> = max.into_iter().map(|max| (max as i64) * offset | 0x1FFFF).collect();
        reverser.add_modulo_measured_seeds(min, max, (bound as i64) * offset);
    }
}

pub fn add_next_int_call(reverser: &mut RandomReverser, min: Vec<i32>, max: Vec<i32>) {
    let offset = 1_i64 << 16;
    let min: Vec<i64> = min.into_iter().map(|min| (min as i64) * offset).collect();
    let max: Vec<i64> = max.into_iter().map(|max| (max as i64) * offset + offset - 1).collect();
    reverser.add_measured_seed(min, max);
}

pub fn add_next_boolean_call(reverser: &mut RandomReverser, value: bool) {
    if value {
        add_bound_next_int_call(reverser, 2, vec![1], vec![1]);
    } else {
        add_bound_next_int_call(reverser, 2, vec![0], vec![0]);
    }
}

pub fn add_next_float_call(reverser: &mut RandomReverser, min: Vec<f32>, max: Vec<f32>) {
    let min_seeds = min
        .into_iter()
        .map(|min| {
            let min_long = (min * (1 << 24) as f32).ceil() as i64;
            min_long << 24
        })
        .collect();

    let max_seeds = max
        .into_iter()
        .map(|max| {
            let max_long = (max * (1 << 24) as f32).ceil() as i64 - 1;
            (max_long << 24) | 0xFFFFFF
        })
        .collect();

    reverser.add_measured_seed(min_seeds, max_seeds);
}

pub fn add_next_long_call(reverser: &mut RandomReverser, min: Vec<i64>, max: Vec<i64>) {
    let min_first_seeds = min
        .iter()
        .map(|min| {
            let min_sign_bit = (min & 0x8000_0000) != 0;
            (min.ushr(32) + if min_sign_bit { 1 } else { 0 }) << 16
        })
        .collect();

    let max_first_seeds = max
        .iter()
        .map(|max| {
            let max_sign_bit = (max & 0x8000_0000) != 0;
            (max.ushr(32) + if max_sign_bit { 2 } else { 1 }) << 16
        })
        .collect();

    reverser.add_measured_seed(min_first_seeds, max_first_seeds);

    let mut second_seeds_min: Vec<i64> = Vec::new();
    let mut second_seeds_max: Vec<i64> = Vec::new();

    for (min, max) in min.into_iter().zip(max.into_iter()) {
        if max - min < 1 << 32 && 0 <= max - min {
            second_seeds_min.push((min & 0xFFFF_FFFF) << 16);
            second_seeds_max.push((((max & 0xFFFF_FFFF) + 1) << 16) - 1);
        }
    }

    reverser.add_measured_seed(second_seeds_min, second_seeds_max);
}

pub fn add_next_double_call(reverser: &mut RandomReverser, min: Vec<f64>, max: Vec<f64>) {
    let min: Vec<i64> = min.into_iter().map(|min| (min * (1_i64 << 53) as f64).ceil() as i64).collect();
    let max: Vec<i64> = max.into_iter().map(|max| (max * (1_i64 << 53) as f64).ceil() as i64 - 1).collect();

    let min_first_seeds = min.iter().map(|min_long| (min_long << 27) << 22).collect();
    let max_first_seeds = max.iter().map(|max_long| ((max_long << 27) << 22) | 0x3FFFFF_i64).collect();

    reverser.add_measured_seed(min_first_seeds, max_first_seeds);

    let mut second_seeds_min: Vec<i64> = Vec::new();
    let mut second_seeds_max: Vec<i64> = Vec::new();

    for (min, max) in min.into_iter().zip(max.into_iter()) {
        if ((max - min) % (1 << 53)) < (1 << 27) {
            second_seeds_min.push((min & 0x7FFFFFF) << 21);
            second_seeds_max.push(((max & 0x7FFFFFF) << 21) | 0x1FFFFF);
        }
    }

    reverser.add_measured_seed(second_seeds_min, second_seeds_max);
}
