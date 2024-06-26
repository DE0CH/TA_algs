use ordered_float::NotNan;
use ta_algs::common::PointsGrid;
use ta_algs::entrance::get_raw_points;
use itertools::Itertools;

fn temperature(current: usize, max: usize) -> NotNan<f64> {
    let current = current as f64;
    let max = max as f64;
    let x = (current + 1.0) / max;
    let a = 7.0 as f64;
    let b = 1 as f64;
    let t = (-a*x).exp() - (-a).exp() + (1.0 - x) * (-a).exp();
    let t = t * b;
    // \left(e^{-ax}-e^{-a}\right)+\left(1-x\right)\cdot e^{-a}
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
fn main() {
    let (points, seed, iterations) = get_raw_points();
    let points_grid = PointsGrid::new(points);
    let mut rng = fastrand::Rng::with_seed(seed);
    let xc = points_grid.generate_random_point(&mut rng);
    let mut xc = xc.round_down();
    let xc_sn = xc.clone();

    let mut global = xc_sn.get_bardelta();
    let mut current = xc_sn.get_bardelta();

    let total_iterations = iterations as usize;
    for current_iteration in 0..total_iterations {
        let t = temperature(current_iteration, total_iterations);
        let k: Vec<_> = points_grid.get_k(current_iteration, total_iterations).collect();
        let mc = points_grid.get_mc(current_iteration, total_iterations);
        let y = xc.generate_neighbor_point(&k, mc, &mut rng);
        let ym = y.round_down();
        let ym_sn = ym.clone();
        let new_current = ym_sn.get_bardelta();
        global = global.max(new_current);
        if pp(current, new_current, t) >= NotNan::new(rng.f64()).unwrap() {
            xc = ym;
        }
        current = new_current;
    }
    println!("{}", global)
}
