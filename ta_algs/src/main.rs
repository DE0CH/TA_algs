pub mod common;
pub mod ta_shirink_delta;
pub mod entrance;
pub mod utils;

fn main() {
    let (points, seed, iterations) = entrance::get_raw_points();
    let mut rng = fastrand::Rng::with_seed(seed);
    ta_shirink_delta::trial(points, iterations, &mut rng);
}

