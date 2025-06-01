use crate::util::Random;

fn next(random: &mut Random, bits: i32) -> i32 {
    random.seed = random.lcg.next_seed(random.seed);
    ((random.seed as u64) >> (48 - bits)) as i32
}

/// An implementation of Java's nextBoolean().
pub fn next_boolean(random: &mut Random) -> bool {
    next(random, 1) == 1
}

/// An implementation of Java's nextInt().
pub fn next_int(random: &mut Random) -> i32 {
    next(random, 32)
}

/// An implementation of Java's nextInt(bound).
/// Will panic if bound is < 1.
pub fn next_bound_int(random: &mut Random, bound: i32) -> i32 {
    if bound < 1 {
        panic!("bound must be positive");
    }

    if (bound & -bound) == bound {
        // Bound is power of 2
        let multiplier: i64 = next(random, 31) as i64;
        return ((bound as i64 * multiplier) >> 31) as i32;
    }

    let mut bits: i32;
    let mut value: i32;

    loop {
        bits = next(random, 31);
        value = bits % bound;

        if bits - value + (bound + 1) >= 0 {
            break;
        }
    }

    value
}

/// An implementation of Java's nextFloat().
pub fn next_float(random: &mut Random) -> f32 {
    next(random, 24) as f32 / ((1 << 24) as f32)
}

/// An implementation of Java's nextLong().
pub fn next_long(random: &mut Random) -> i64 {
    let upper_bits = (next(random, 32) as i64) << 32;
    upper_bits + next(random, 32) as i64
}

/// An implementation of Java's nextDouble().
pub fn next_double(random: &mut Random) -> f64 {
    let upper_bits = (next(random, 26) << 27) as i64;
    (upper_bits + next(random, 27) as i64) as f64 / (1_i64 << 53) as f64
}