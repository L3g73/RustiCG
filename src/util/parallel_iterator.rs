use std::sync::mpsc;
use std::sync::mpsc::{Receiver, Sender};

/// An iterator whose sources can be in multiple threads.
pub struct ParallelIterator<T> {
    sender: Option<Sender<T>>,
    receiver: Receiver<T>,
}

impl<T> ParallelIterator<T> {
    pub(crate) fn new() -> ParallelIterator<T> {
        let (sender, receiver) = mpsc::channel();
        ParallelIterator {
            sender: Some(sender),
            receiver,
        }
    }

    pub(crate) fn get_new_sender(&self) -> Sender<T> {
        if let Some(sender) = &self.sender {
            return sender.clone();
        }

        panic!("Called get_new_sender() after the iterator has started")
    }
}

impl<T> Iterator for ParallelIterator<T> {
    type Item = T;

    fn next(&mut self) -> Option<Self::Item> {
        if self.sender.is_some() {
            self.sender = None; // Drop reference
        }

        self.receiver.recv().ok()
    }
}
