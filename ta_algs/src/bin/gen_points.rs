use fastrand::Rng;
use std::iter::repeat_with;
use itertools::Itertools;
use clap::{arg, value_parser, Command};

fn choose_from_array<'a, T>(v: &'a Vec<T>, rng: &mut Rng) -> &'a T {
    let i = rng.usize(0..v.len());
    &v[i]
}

fn main() {
    let matches = Command::new("")
        .arg(
            arg!(dimension: <dimension> "Dimension")
            .value_parser(value_parser!(usize))
        ).arg(
            arg!(n_points: <n_points> "Number of points")
            .value_parser(value_parser!(usize))
        ).arg(
            arg!(seed: -s --seed [seed] "Seed for the rng")
            .value_parser(value_parser!(u64))
            .default_value("42")
        ).get_matches();
    let seed = *matches.get_one::<u64>("seed").unwrap();
    let d = *matches.get_one::<usize>("dimension").unwrap();
    let n_points = *matches.get_one::<usize>("n_points").unwrap();
    let mut rng = Rng::with_seed(seed);
    let points: Vec<_> = repeat_with(|| {
        let t: Vec<i64> = repeat_with(|| rng.i64(0..1000)).take(n_points).collect();
        t
    }).take(d).collect();
    for i in 0..n_points {
        let s = (0..d).map(|id| points[id][i] as f64 / 1000.0).map(|x| format!("{:.3}", x));
        println!("{}", s.format(" "));
    }
}
