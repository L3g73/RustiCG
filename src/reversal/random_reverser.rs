use crate::math::component::big_matrix::BigMatrix;
use crate::math::component::big_vector::{OwnedBigVector, Vector};
use crate::math::lattice::decomposition::inverse;
use crate::math::lattice::enumerate::enumerate;
use crate::math::lattice::lll;
use crate::util::util::LCM;
use crate::util::util::{filled_vec, ToU64};
use crate::util::{ParallelIterator, Random, JAVA, LCG};
use malachite::base::num::arithmetic::traits::Pow;
use malachite::base::num::basic::traits::{One, Two, Zero};
use malachite::base::num::conversion::traits::ToSci;
use malachite::rational::Rational;
use malachite::Integer;

/// Skips [LCG] calls and applies a filter function on the resulting seeds instead.<br>
/// Is usually not interacted with directly, use [DynamicProgram] instead.
pub struct FilteredSkip {
    skip_lcg: LCG,
    filter: Box<dyn Fn(&mut Random) -> bool>,
}

impl FilteredSkip {
    /// Creates a new FilteredSkip.
    pub fn new(current_index: i64, filter: Box<dyn Fn(&mut Random) -> bool>) -> Self {
        FilteredSkip {
            skip_lcg: JAVA.combine(current_index),
            filter,
        }
    }

    /// Checks whether this filter matches the `rand`.
    pub fn check_state(&self, rand: &mut Random) -> bool {
        rand.advance(&self.skip_lcg);
        (self.filter)(rand)
    }
}

#[derive(Clone)]
struct Bounds {
    mins: Vec<Integer>,
    maxes: Vec<Integer>,
}

/// Generates the lattice used and applies bounds / filters.<br>
/// Is usually not interacted with directly, use [DynamicProgram] instead.<br>
/// See [here](https://gist.github.com/EDDxample/38a9acddcd29f15af034fd91da93b8fa) for an explanation on how the lattice works.
///
/// # Bounds vs Filters
/// Bounds add a bound to the seed region in the lattice.<br>
/// Filters are a more brute-force way of reducing seeds, being a predicate function called for each resulting seed.
/// ## Which one should I use?
/// For strong restrictions like `number % 16 == 0` [^example_fn] a bound is usually the best way.<br>
/// Inverted restrictions like `number % 16 != 0` will be simplified to `number % 16 > 0`, they can be used as bounds without problems as well.<br>
/// Restrictions that result in two bounds[^two_bound_restriction] will also result in multiple lattices having to be generated, causing a big performance impact.<br>
/// For those cases, it's recommended to use a filter instead.
/// <br>
/// There's no general answer, as often restrictions can be more efficient if they're filters, even when they only result in one lattice.<br>
/// It's recommended to prefer stronger restrictions[^strong_restriction] as bounds over weaker ones[^weak_restriction].
/// But too few bounds will also result in a worse performance, so **it has to be decided on a case by case basis**.
/// 
/// [^example_fn]: All examples use the [next_bound_int()](java::next_bound_int) fn.
/// [^two_bound_restriction]: e.g. `number % 16 != 5` used as a bound will internally be represented as `number % 16 <= 4 || number % 16 >= 6`
/// [^strong_restriction]: e.g. `number % 16 == 0`, will remove 93,75% of all seeds
/// [^weak_restriction]: e.g. `number % 16 > 3`, will only remove 25%
pub struct RandomReverser {
    modulus: Integer,
    multiplier: Integer,

    lcg: LCG,
    bounds: Vec<Bounds>,
    call_indices: Vec<i64>,
    filtered_skips: Vec<FilteredSkip>,
    lattice: BigMatrix,
    current_call_index: i64,
    dimensions: usize,
    success_chance: Rational
}

impl RandomReverser {
    /// Creates a new RandomReverser with the given LCG.
    pub fn new(lcg: LCG) -> RandomReverser {
        let mut modulus = Integer::from(lcg.modulus);

        if lcg.modulus <= 0 {
            modulus += Integer::TWO.pow(64);
        }

        RandomReverser {
            multiplier: Integer::from(lcg.multiplier) % &modulus,
            modulus,
            lcg,
            bounds: vec![Bounds {
                mins: Vec::new(),
                maxes: Vec::new(),
            }],
            call_indices: Vec::new(),
            filtered_skips: Vec::new(),
            lattice: BigMatrix::new(0, 0),
            current_call_index: 0,
            dimensions: 0,
            success_chance: Rational::ONE
        }
    }

    /// Starts reversing seeds.
    pub fn into_all_valid_seeds(self) -> Box<dyn Iterator<Item = i64>> {
        if self.success_chance != Rational::ONE {
            eprintln!("Ignored ~{}% of all seeds", (Rational::ONE - &self.success_chance).to_sci());
        }

        let dim = self.dimensions;

        let seeds: Box<dyn Iterator<Item = i64>> = if dim == 0 {
            Box::new(0..self.lcg.modulus)
        } else {
            let lattices = self.get_reduced_lattices();
            let lower = filled_vec(self.bounds.len(), |i| {
                let bounds = &self.bounds[i];
                OwnedBigVector::new_from(dim, |i| Rational::from(&bounds.mins[i]))
            });
            let upper = filled_vec(self.bounds.len(), |i| {
                let bounds = &self.bounds[i];
                OwnedBigVector::new_from(dim, |i| Rational::from(&bounds.maxes[i]))
            });
            let mut offset = OwnedBigVector::new(dim);
            let mut rand = Random::of_seed(self.lcg.clone(), 0);

            for i in 0..dim {
                offset.set(
                    i,
                    Rational::from_integers(rand.get_seed().into(), Integer::ONE),
                );

                if i != dim - 1 {
                    rand.advance_steps(self.call_indices[i + 1] - self.call_indices[i]);
                }
            }

            let seed_reverser = self.lcg.combine(-self.call_indices[0]);

            let iterator = ParallelIterator::new();

            for ((lattice, lower), upper) in lattices
                .into_iter()
                .zip(lower.into_iter())
                .zip(upper.into_iter())
            {
                let basis = lattice.transpose();
                let sender = iterator.get_new_sender();
                enumerate(lattice.transpose(), lower, upper, &offset, move |v| {
                    sender.send(&basis * v).is_err()
                });
            }

            Box::new(
                iterator
                    .map(move |vec| vec + &offset)
                    .map(|v| v.get(0).numerator_ref().to_u64() as i64)
                    .map(move |l| seed_reverser.next_seed(l)),
            )
        };

        Box::new(seeds.filter(move |seed| {
            for call in &self.filtered_skips {
                let mut rr: Random = Random::of_seed(self.lcg.clone(), *seed);
                if !call.check_state(&mut rr) {
                    return false;
                }
            }

            true
        }))
    }

    fn get_reduced_lattices(&self) -> Vec<BigMatrix> {
        let dim = self.dimensions;
        let mut lattices = Vec::with_capacity(self.bounds.len());

        for bounds in &self.bounds {
            let mut side_length: Vec<Integer> = Vec::with_capacity(dim);
            for i in 0..dim {
                side_length.push(&bounds.maxes[i] - &bounds.mins[i] + Integer::ONE);
            }

            let mut lcm_value = Integer::ONE;
            for i in 0..dim {
                lcm_value = (&lcm_value).lcm(&side_length[i]);
            }

            let mut scales = BigMatrix::new(dim, dim);
            for i in 0..dim {
                scales.set(i, i, Rational::from(&lcm_value / &side_length[i]));
            }

            let scaled_lattice = &self.lattice * &scales;
            let reduced = lll::reduce(scaled_lattice);
            lattices.push(reduced * inverse(&scales));
        }

        lattices
    }

    fn add_dimensions(&mut self, added_dimensions: usize) {
        self.dimensions += added_dimensions;
        let mut new_lattice = BigMatrix::new(self.dimensions + 1, self.dimensions);

        // Copy old lattice
        if self.lattice.column_count != 0 {
            for row in 0..self.dimensions + 1 - added_dimensions {
                for col in 0..self.dimensions - added_dimensions {
                    new_lattice.set(row, col, self.lattice.get(row, col).clone())
                }
            }
        }

        self.lattice = new_lattice;
    }

    /// Adds a [FilteredSkip].
    pub fn add_filtered_skip(&mut self, filtered_skip: FilteredSkip) {
        self.filtered_skips.push(filtered_skip);
    }

    fn add_bounds(&mut self, min: Vec<Integer>, max: Vec<Integer>) {
        let original_len = self.bounds.len();
        for _ in 0..min.len() - 1 {
            self.bounds.extend_from_within(0..original_len);
        }

        for (i, bounds) in self.bounds.iter_mut().enumerate() {
            bounds.mins.push(min[i / original_len].clone());
            bounds.maxes.push(max[i / original_len].clone());
        }
    }

    /// Skips the given amount of [LCG] calls.
    pub fn add_unmeasured_seeds(&mut self, num_seeds: i64) {
        self.current_call_index += num_seeds;
    }

    /// Adds bounds to the next seed.<br>
    /// `mins` and `maxes` are the minimum and maximum bounds of the ranges.<br>
    /// The resulting seeds will be in at least one of the ranges.<br>
    /// Multiple ranges are not recommended unless absolutely necessary due to performance reasons, see [RandomReverser] for more information.
    pub fn add_measured_seed<T: Into<Integer>>(&mut self, mins: Vec<T>, maxes: Vec<T>) {
        if mins.len() != maxes.len() {
            panic!("The amount of min seeds is not equal to the amount of max seeds!");
        }

        if mins.len() == 0 {
            self.add_unmeasured_seeds(1);
        }

        let min: Vec<Integer> = mins.into_iter().map(|m| m.into() % &self.modulus).collect();

        let max: Vec<Integer> = maxes
            .into_iter()
            .enumerate()
            .map(|(i, m)| {
                let max = m.into() % &self.modulus;
                if max < min[i] {
                    max + &self.modulus
                } else {
                    max
                }
            })
            .collect();

        self.add_bounds(min, max);

        self.current_call_index += 1;
        self.call_indices.push(self.current_call_index);

        self.add_dimensions(1);
        let dims = self.dimensions;

        let mult = (&self.multiplier)
            .pow((self.call_indices[dims - 1] - self.call_indices[0]) as u64)
            % &self.modulus;
        self.lattice.set(0, dims - 1, Rational::from(mult));
        self.lattice
            .set(dims, dims - 1, Rational::from(self.modulus.clone()));
    }

    /// Adds bounds to the next seed.<br>
    /// `mins` and `maxes` are the minimum and maximum bounds of the seed ranges, modulo `measured_mod`.<br>
    /// See [RandomReverser::add_measured_seed]
    pub fn add_modulo_measured_seeds<T: Into<Integer>>(
        &mut self,
        min: Vec<T>,
        max: Vec<T>,
        measured_mod: T,
    ) {
        let measured_mod = measured_mod.into();

        let min: Vec<Integer> = min
            .into_iter()
            .map(|min| min.into() % &measured_mod)
            .collect();

        let max: Vec<Integer> = max
            .into_iter()
            .enumerate()
            .map(|(i, max)| {
                let max = max.into() % &measured_mod;
                if max < min[i] {
                    max + &measured_mod
                } else {
                    max
                }
            })
            .collect();

        let residue = &self.modulus % &measured_mod;
        if residue == 0 {
            self.add_bounds(min, max);
            self.current_call_index += 1;
            self.call_indices.push(self.current_call_index);

            self.add_dimensions(1);

            let dims = self.dimensions;
            if dims == 1 && self.modulus != measured_mod {
                eprintln!("First call not a bound on a seed. Junk output may be produced.");
            }

            let mult = (&self.multiplier)
                .pow((self.call_indices[dims - 1] - self.call_indices[0]) as u64)
                % &self.modulus;
            self.lattice.set(0, dims - 1, Rational::from(mult));
            self.lattice
                .set(dims, dims - 1, Rational::from(measured_mod));
            return;
        }

        self.success_chance *= Rational::ONE - Rational::from_integers_ref(&residue, &self.modulus);

        self.add_bounds(vec![Integer::ZERO], vec![&self.modulus - residue]); // seeds not fitting in this range are currently unsupported
        self.current_call_index += 1;
        self.call_indices.push(self.current_call_index);

        self.add_bounds(min, max);
        self.call_indices.push(self.current_call_index);

        self.add_dimensions(2);

        let dims = self.dimensions;
        let mult = (&self.multiplier)
            .pow((self.call_indices[dims - 1] - self.call_indices[0]) as u64)
            % &self.modulus;
        self.lattice.set(0, dims - 2, Rational::from(&mult));
        self.lattice.set(0, dims - 1, Rational::from(&mult));

        // Vector capturing the effect of the modulo 2^48 operation on the residue class modulo measured_mod
        self.lattice
            .set(dims - 1, dims - 1, Rational::from(&self.modulus));
        self.lattice
            .set(dims - 1, dims - 2, Rational::from(&self.modulus));

        // Vector identifying everything in residue classes modulo measured_mod
        self.lattice
            .set(dims, dims - 1, Rational::from(measured_mod))
    }
}
