pub mod common;
pub mod ta_shrink_delta;
pub mod ta_shrink_deltabar;
pub mod entrance;
pub mod utils;

use common::PointsGrid;
use common::get_index_up;
use common::get_index_down;
use common::Index;

fn main() {
    let (points, seed, iterations) = entrance::get_raw_points();
    let mut rng = fastrand::Rng::with_seed(seed);
    let points2 = points.clone();
    let points = PointsGrid::new(points);
    let (a, b, c) = ta_shrink_deltabar::trial(points, iterations, &mut rng);
    println!("The best point is {:?} with a global delta of {:?} and a current delta of {:?}", a, b, c);
    let coordinate = a.into_iter().enumerate().map(|(id, i)| {
        match i {
            Index::Zero => 0.0,
            Index::One => 1.0,
            Index::N(x) => points2[id][x].into(),
        }
    }).collect::<Vec<_>>();
    println!("The best point in the original coordinate system is {:?}", coordinate);
}

