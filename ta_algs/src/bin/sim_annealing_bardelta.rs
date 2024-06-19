use ordered_float::NotNan;
use ta_algs::common::PointsGrid;
use ta_algs::entrance::get_raw_points;

fn temperature(current: usize, max: usize) -> NotNan<f64> {
    let current = current as f64;
    let max = max as f64;
    let t = 1.0 - (current + 1.0) / max;
    NotNan::new(t).unwrap()
}

fn pp(current: NotNan<f64>, new: NotNan<f64>, t: NotNan<f64>) -> NotNan::<f64> {
    if new > current {
        return NotNan::new(1.0).unwrap()
    }
    if t == 0.0 {
        return NotNan::new(0.0).unwrap()
    }

    let p = (new - current) / t;
    NotNan::new(p.exp()).unwrap()
}
// Let s = s0
// For k = 0 through kmax (exclusive):
// T ← temperature( 1 - (k+1)/kmax )
// Pick a random neighbour, snew ← neighbour(s)
// If P(E(s), E(snew), T) ≥ random(0, 1):
// s ← snew
// Output: the final state s
fn main() {
    let (points, seed, iterations) = get_raw_points();
    let points_grid = PointsGrid::new(points);
    let mut rng = fastrand::Rng::with_seed(seed);
    let xc = points_grid.generate_random_point(&mut rng);
    let mut xc = xc.round_down();

    let mut global = xc.get_bardelta();

    let total_iterations = iterations as usize;
    for current_iteration in 0..total_iterations {
        let t = temperature(current_iteration, total_iterations);
        let k: Vec<_> = points_grid.get_k(current_iteration, total_iterations).collect();
        let mc = points_grid.get_mc(current_iteration, total_iterations);
        let y = xc.generate_neighbor_point(&k, mc, &mut rng);
        let ym = y.round_down();
        global = global.max(ym.get_bardelta());

        if pp(xc.get_bardelta(), ym.get_bardelta(), t) >= NotNan::new(rng.f64()).unwrap() {
            xc = ym;
        }
    }
    println!("{}", global)
}
