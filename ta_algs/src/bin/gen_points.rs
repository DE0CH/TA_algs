use fastrand::Rng;
use std::iter::repeat_with;
use itertools::Itertools;

fn choose_from_array<'a, T>(v: &'a Vec<T>, rng: &mut Rng) -> &'a T {
    let i = rng.usize(0..v.len());
    &v[i]
}

fn main() {
    let mut rng = Rng::new();
    let n_points = 30;
    let d = 3;
    let choices: Vec<_> = repeat_with(|| rng.i64(0..1000)).take(20).collect();
    let points: Vec<_> = repeat_with(|| {
        let t: Vec<i64> = repeat_with(|| choose_from_array(&choices, &mut rng)).take(n_points).copied().collect();
        t
    }).take(d).collect();
    for i in 0..n_points {
        let s = (0..d).map(|id| points[id][i] as f64 / 1000.0).map(|x| format!("{:.3}", x));
        println!("{}", s.format(" "));
    }
}
