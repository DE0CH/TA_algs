use ta_algs::entrance::get_raw_points;
use ta_algs::common::PointsGrid;
use fastrand::Rng;
use itertools::Itertools;


fn main() {
    let mut rng = Rng::with_seed(42);
    let (raw_points, _, _) = get_raw_points();
    let points_grid = PointsGrid::new(raw_points);
    let point = points_grid.generate_random_point(&mut rng);
    println!("random point: {}", point.float().map(|x| format!("{:.3}", *x)).format(" "));
    let point = point.round_down();
    println!("rounded down: {}", point.float().map(|x| format!("{:.3}", *x)).format(" "));
    let point = point.snap_down();
    println!("snapped down: {}", point.float().map(|x| format!("{:.3}", *x)).format(" "));
}
