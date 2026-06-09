use crate::math::component::big_matrix::BigMatrix;
use crate::math::component::big_vector::{OwnedBigVector, Vector, ViewBigVector};
use crate::math::lattice::decomposition::inverse;
use crate::math::optimize::optimize::{Optimize, OptimizeBuilder};
use malachite::Integer;
use malachite::base::num::arithmetic::traits::{Ceiling, Floor};
use malachite::base::num::basic::traits::One;
use malachite::rational::Rational;
use std::sync::atomic::{AtomicUsize, Ordering};
use std::sync::{Arc, Mutex};
use std::thread::JoinHandle;

pub fn enumerate<Consumer: Fn(OwnedBigVector) -> bool + Sync + Send + 'static>(
    basis: BigMatrix,
    lower: OwnedBigVector,
    upper: OwnedBigVector,
    origin: &OwnedBigVector,
    consumer: Consumer,
) {
    let mut builder = OptimizeBuilder::of_size(basis.row_count);

    let mut lower_into = lower.numbers.into_iter();
    let mut upper_into = upper.numbers.into_iter();

    for i in 0..basis.row_count {
        builder.with_lower_bound(i, lower_into.next().unwrap());
        builder.with_upper_bound(i, upper_into.next().unwrap());
    }

    let constraints = builder.build();
    let root_inverse = inverse(&basis);
    let root_origin = &root_inverse * origin;

    let root_size = basis.row_count;
    let root_fixed = OwnedBigVector::new(root_size);
    let root_constraints = constraints.clone();

    let mut widths: Vec<Rational> = Vec::new();
    let mut order: Vec<usize> = Vec::new();

    for i in 0..root_size {
        let gradient = &root_inverse.get_row(i);
        let min = constraints.clone().minimize(gradient);
        let max = constraints.clone().maximize(gradient);

        widths.push(max - min);
        order.push(i);
    }

    order.sort_by_key(|i| &widths[*i]);

    let root = SearchNode {
        size: root_size,
        depth: 0,
        fixed: root_fixed,
        constraints: root_constraints,
    };

    root.compute(root_inverse, root_origin, order, consumer);
}

struct SearchNode {
    pub size: usize,
    pub depth: usize,
    pub fixed: OwnedBigVector,
    pub constraints: Optimize,
}

impl SearchNode {
    fn create_child(
        &self,
        index: usize,
        i: Integer,
        gradient: ViewBigVector,
        offset: &Rational,
    ) -> SearchNode {
        let i = Rational::from(i);
        let next_optimize = self.constraints.with_strict_bound(&gradient, &i + offset);
        let mut next_fixed = self.fixed.clone();
        let updated = next_fixed.get(index) + &i;
        next_fixed.set(index, updated);

        SearchNode {
            size: self.size,
            depth: self.depth + 1,
            fixed: next_fixed,
            constraints: next_optimize,
        }
    }

    pub fn compute<Consumer: Fn(OwnedBigVector) -> bool + Sync + Send + 'static>(
        self,
        inverse: BigMatrix,
        origin: OwnedBigVector,
        order: Vec<usize>,
        consumer: Consumer,
    ) {
        let threads: usize = std::thread::available_parallelism().unwrap().get();

        let inverse = Arc::new(inverse);
        let origin = Arc::new(origin);
        let order = Arc::new(order);
        let consumer = Arc::new(consumer);
        let stack = Arc::new(Mutex::new(vec![self]));
        let waiting_threads = Arc::new(AtomicUsize::new(threads));

        let mut handles = vec![];
        for _ in 0..threads {
            handles.push(Self::spawn_compute_thread(
                threads,
                stack.clone(),
                consumer.clone(),
                waiting_threads.clone(),
                inverse.clone(),
                origin.clone(),
                order.clone(),
            ))
        }
    }

    fn spawn_compute_thread<Consumer: Fn(OwnedBigVector) -> bool + Sync + Send + 'static>(
        threads: usize,
        stack: Arc<Mutex<Vec<SearchNode>>>,
        consumer: Arc<Consumer>,
        waiting_threads: Arc<AtomicUsize>,
        inverse: Arc<BigMatrix>,
        origin: Arc<OwnedBigVector>,
        order: Arc<Vec<usize>>,
    ) -> JoinHandle<()> {
        std::thread::spawn(move || {
            'main_loop: loop {
                while let Some(node) = { stack.lock().unwrap().pop() } {
                    waiting_threads.fetch_sub(1, Ordering::SeqCst);
                    if node.depth == node.size {
                        if consumer(node.fixed) {
                            break 'main_loop; // Receiver has been dropped
                        }
                        waiting_threads.fetch_add(1, Ordering::SeqCst);
                        continue;
                    }

                    let index = order[node.depth];

                    let gradient = &inverse.get_row(index);
                    let offset = origin.get(index);

                    let mut min = (node.constraints.clone().minimize(gradient) - offset).ceiling();
                    let max = (node.constraints.clone().maximize(gradient) - offset).floor();

                    while min <= max {
                        let child = node.create_child(
                            index,
                            min.clone(),
                            inverse.get_row(index),
                            origin.get(index),
                        );
                        {
                            stack.lock().unwrap().push(child);
                        }
                        min += Integer::ONE;
                    }
                    waiting_threads.fetch_add(1, Ordering::SeqCst);
                }

                if waiting_threads.load(Ordering::Relaxed) == threads {
                    break;
                }
            }
        })
    }
}
