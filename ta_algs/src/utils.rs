use std::mem::ManuallyDrop;
use std::mem::size_of;
use std::ptr;
use std::alloc;
use std::alloc::Layout;
use itertools::Itertools;
use itertools::Tee;

pub struct PermutationIterator<'a, T> {
    vec: ManuallyDrop<Vec<T>>,
    start: usize,
    end: usize,
    permutation: &'a [usize],
}

impl<'a, T> Iterator for PermutationIterator<'a, T> {
    type Item = T;

    fn next(&mut self) -> Option<Self::Item> {
        if self.start == self.end {
            None
        } else {
            let j = self.permutation[self.start];
            let p = &self.vec[j] as *const T;
            let result = unsafe { ptr::read(p) };
            self.start += 1;
            Some(result)
        }
    }

    fn size_hint(&self) -> (usize, Option<usize>) {
        let len = (self.end as usize - self.start as usize)
                  / size_of::<T>();
        (len, Some(len))
    }
}

impl<'a, T> DoubleEndedIterator for PermutationIterator<'a, T> {
    fn next_back(&mut self) -> Option<T> {
        if self.start == self.end {
            None
        } else {
            let j = self.permutation[self.end - 1];
            let p = &self.vec[j] as *const T;
            let result = unsafe { ptr::read(p) };
            self.end -= 1;
            Some(result)
        }
    }
}

impl<'a, T> Drop for PermutationIterator<'a, T> {
    fn drop(&mut self) {
        if self.vec.capacity() == 0 {
            return;
        }
        for _ in &mut *self {}
        let layout = Layout::array::<T>(self.vec.capacity()).unwrap();
        unsafe {
            alloc::dealloc(self.vec.as_ptr() as *mut u8, layout);
        }
    }
}

pub fn into_iterator_with_permutation<'a, T: 'a>(
    stuff: Vec<T>,
    permutation: &'a [usize],
) -> PermutationIterator<'a, T> {
    assert_eq!(stuff.len(), permutation.len(), "Permutation length does not match the length of the vector");
    let mut check_array = vec![false; stuff.len()];
    permutation.iter().for_each(|x| {
        assert!(!check_array[*x], "Permutation is not a permutation");
        check_array[*x] = true;
    });
    unsafe { into_iterator_with_permutation_unchecked(stuff, permutation) }
}

pub unsafe fn into_iterator_with_permutation_unchecked<'a, T: 'a>(
    stuff: Vec<T>,
    permutation: &'a [usize],
) -> PermutationIterator<'a, T> {
    let len = stuff.len();
    let vec = PermutationIterator {
        vec: ManuallyDrop::new(stuff),
        start: 0,
        end: len,
        permutation,
    };
    vec
}

pub fn print_iterator<I: Iterator<Item = impl std::fmt::Debug + Clone>>(msg: &str, iterator: I) -> Tee<I>{
    let (a, b) = iterator.tee();
    println!("{} {:?}", msg, b.collect::<Vec<_>>());
    a
}
