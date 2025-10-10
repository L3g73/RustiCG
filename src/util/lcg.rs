use crate::util::util::Ushr;

/// [The LCG used by Java in versions prior to 17.](https://hg.openjdk.org/jdk8/jdk8/jdk/file/tip/src/share/classes/java/util/Random.java#l88)
pub const JAVA: LCG = LCG::new(0x5DEECE66D, 0xB, 1 << 48);

/// A [linear conguential generator](https://en.wikipedia.org/wiki/Linear_congruential_generator).
#[derive(Eq, PartialEq, Clone, Debug)]
pub struct LCG {
    pub multiplier: i64,
    pub addend: i64,
    pub modulus: i64,
    can_mask: bool,
}

impl LCG {

    /// Creates a new LCG.
    pub const fn new(multiplier: i64, addend: i64, modulus: i64) -> LCG {
        LCG {
            multiplier,
            addend,
            modulus,
            can_mask: (modulus & -modulus) == modulus,
        }
    }

    pub fn next_seed(&self, seed: i64) -> i64 {
        self.modulus(
            seed.overflowing_mul(self.multiplier).0
                .overflowing_add(self.addend).0,
        )
    }

    fn modulus(&self, n: i64) -> i64 {
        if self.can_mask {
            n & (self.modulus - 1)
        } else {
            let n_unsiged = u64::from_ne_bytes(n.to_ne_bytes());
            let mod_unsigned = u64::from_ne_bytes(self.modulus.to_ne_bytes());
            i64::from_ne_bytes((n_unsiged % mod_unsigned).to_ne_bytes())
        }
    }

    pub fn combine(&self, steps: i64) -> LCG {
        let mut multiplier: i64 = 1;
        let mut addend: i64 = 0;

        let mut intermediate_multiplier = self.multiplier;
        let mut intermediate_addend = self.addend;

        let mut k = steps;
        while k != 0 {
            if k & 1 != 0 {
                multiplier = multiplier.wrapping_mul(intermediate_multiplier);
                addend = intermediate_multiplier
                    .wrapping_mul(addend)
                    .wrapping_add(intermediate_addend);
            }

            intermediate_addend = (intermediate_multiplier + 1).wrapping_mul(intermediate_addend);
            intermediate_multiplier = intermediate_multiplier.wrapping_mul(intermediate_multiplier);

            k = k.ushr(1);
        }

        multiplier = self.modulus(multiplier);
        addend = self.modulus(addend);

        LCG::new(multiplier, addend, self.modulus)
    }
}

/// A generator for random numbers.
pub struct Random {
    /// The [LCG] used.
    pub lcg: LCG,

    /// The current internal seed.
    pub seed: i64
}

impl Random {

    /// Creates a new random generator.
    pub fn of_seed(lcg: LCG, seed: i64) -> Random {
        Random { lcg, seed }
    }

    /// Returns the internal seed.
    pub fn get_seed(&self) -> i64 {
        self.seed
    }

    /// Sets the internal seed.
    pub fn set_seed(&mut self, seed: i64) {
        self.seed = seed;
    }

    /// Advances the internal seed by the given number of steps.
    pub fn advance_steps(&mut self, steps: i64) {
        self.advance(&self.lcg.combine(steps))
    }

    /// Advances the internal seed by one step using the given LCG.
    pub fn advance(&mut self, skip: &LCG) {
        self.seed = skip.next_seed(self.seed);
    }

}
