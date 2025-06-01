use crate::reversal::random_reverser::{FilteredSkip, RandomReverser};
use crate::util::{Random, LCG};
use std::fmt::Debug;
use std::fmt::Display;
use std::fmt::Formatter;

/// A high-level wrapper for adding lattice bounds.
pub trait CallType {
    fn num_calls(&self) -> i64;

    fn add_to(&self, reverser: &mut RandomReverser) -> Result<(), ReversalError>;
}

/// Converts `self` to a callable function used for filtering the resulting seeds.
pub trait AsFilterFn {
    fn num_calls(&self) -> i64;

    fn as_filter_fn(&self) -> Result<Box<dyn Fn(&mut Random) -> bool>, ReversalError>;
}

/// Allows inverting `self`.
pub trait Invertible {
    fn invert(self) -> Box<Self>;
}

/// A high-level frontend for building the [RandomReverser].
pub struct DynamicProgram {
    random_reverser: RandomReverser,
    current_skip: i64,
    current_index: i64,
}

impl DynamicProgram {

    /// Creates a new DynamicProgram with the given LCG.
    pub fn create(lcg: LCG) -> Self {
        DynamicProgram {
            random_reverser: RandomReverser::new(lcg),
            current_skip: 0,
            current_index: 0,
        }
    }

    /// Skips the given amount LCG calls.
    pub fn skip(&mut self, steps: i64) {
        self.current_skip += steps;
        self.current_index += steps;
    }

    /// Adds a [CallType] to the lattice.<br>
    pub fn add(&mut self, call: Box<dyn CallType>) -> Result<(), ReversalError> {
        self.random_reverser.add_unmeasured_seeds(self.current_skip);
        call.add_to(&mut self.random_reverser)?;
        self.current_skip = 0;
        self.current_index += call.num_calls();
        Ok(())
    }

    /// Adds a [AsFilterFn] to the [RandomReverser].<br>
    pub fn add_filter(&mut self, call: Box<dyn AsFilterFn>) -> Result<(), ReversalError> {
        self.add_boxed_filter_fn(call.as_filter_fn()?, call.num_calls());
        Ok(())
    }

    /// Adds a custom seed filter function to the [RandomReverser].
    pub fn add_filter_fn<F: Fn(&mut Random) -> bool + 'static>(&mut self, filter: F, steps: i64) {
        self.add_boxed_filter_fn(Box::new(filter), steps);
    }

    fn add_boxed_filter_fn(&mut self, filter: Box<dyn Fn(&mut Random) -> bool>, steps: i64) {
        self.random_reverser
            .add_filtered_skip(FilteredSkip::new(self.current_index, filter));
        self.current_skip += steps;
        self.current_index += steps;
    }

    /// Starts reversing seeds.
    pub fn reverse(self) -> impl Iterator<Item = i64> {
        self.random_reverser.into_all_valid_seeds()
    }
}

/// A generic error when processing a [CallType]/[FilterFn](AsFilterFn).
#[derive(Debug)]
pub enum ReversalError {
    RangeIsEmpty,
    TooSmallBound,
    MaxLowerThanMin,
    MinLessThanLowerBound,
    MaxGreaterThanUpperBound,
}

macro_rules! impl_error {
    ($name: ident) => {
        impl Display for $name {
            fn fmt(&self, f: &mut Formatter<'_>) -> std::fmt::Result {
                Debug::fmt(self, f)
            }
        }

        impl std::error::Error for $name {}
    };
}
impl_error!(ReversalError);
