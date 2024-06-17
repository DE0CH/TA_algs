use itertools::Itertools;
use itertools::Tee;

pub fn print_iterator<I: Iterator<Item = impl std::fmt::Debug + Clone>>(msg: &str, iterator: I) -> Tee<I>{
    let (a, b) = iterator.tee();
    println!("{} {:?}", msg, b.collect::<Vec<_>>());
    a
}
