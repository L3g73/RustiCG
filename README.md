# RustiCG
[LattiCG](https://github.com/mjtb49/LattiCG) rewritten in Rust with a focus on speed and usability.

Reverses possible internal seeds of an LCG using a system of inequalities on the output of random calls.<br>
Currently only Java's java.util.Random class is implemented, but other LCGs can be added easily.

Click [here](https://github.com/L3g73/RustiCG/blob/main/CHANGELOG.md) for a changelog.

### Example
```rust 
use rusticg::reversal::dynamic_program::{DynamicProgram, Invertible, ReversalError};
use rusticg::reversal::java;
use rusticg::reversal::java::java_random;
use rusticg::util::{Random, JAVA};

fn main() -> Result<(), ReversalError> {
    let mut program = DynamicProgram::create(JAVA);

    // nextInt() == 0
    program.add(java::next_int().equal_to(0)?)?;

    // nextFloat() <= 0.5
    program.add(java::next_float().less_than_eq(0.5)?)?;

    // nextInt(16) == 2
    program.add(java::next_bound_int(16)?.equal_to(2)?)?;

    // Ignore one call, e.g. an unknown nextInt() call
    program.skip(1);

    // nextInt(4) != 1
    program.add_filter(java::next_bound_int(4)?.equal_to(1)?.invert())?;

    let mut iterator = program.reverse();

    for _ in 0..25 {
        let mut random = Random::of_seed(JAVA, iterator.next().unwrap());
        let r = &mut random;
        print!("{} ", java_random::next_int(r));
        print!("{:.10} ", java_random::next_float(r));
        print!("{} ", java_random::next_bound_int(r, 16));
        java_random::next_int(r); // ignore one call
        println!("{}", java_random::next_bound_int(r, 4));
    }

    Ok(())
}
```