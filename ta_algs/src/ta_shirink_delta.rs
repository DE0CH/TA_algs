use ordered_float::NotNan;
use crate::common::PointsGrid;
use fastrand::Rng;
use itertools::Itertools;

pub fn trial(point_set: Vec<Vec<NotNan<f64>>>, iterations: u64, rng: &mut Rng) {
    let points_grid = PointsGrid::new(point_set);
    let i_tilde = (iterations as f64).sqrt();
    let mut thresh: Vec<_> = (0..i_tilde as usize).map(|i| {
        let current_iteration = (i+1) as f64;
        let k = points_grid.ordered_points.iter().map(|points| {
            let v = (points.len()-1)/2;
            v as f64 * (i_tilde-current_iteration)/i_tilde + current_iteration/i_tilde
        });
        let k: Vec<_> = k.map(|v| v as usize).collect();
        let mc: usize = (2.0 + current_iteration/i_tilde*(points_grid.d as f64 - 2.0)) as usize;
        let xc_index = points_grid.generate_xc_delta(rng);
        let current = xc_index.snap_up(rng).get_delta();
        let xc_plus = xc_index.generate_neighbor_delta(&k, mc, rng);
        let fxc = xc_plus.snap_up(rng).get_delta();
        NotNan::new((fxc-current).abs()).unwrap()
    }).collect();
    thresh.sort();

    let xc = points_grid.generate_xc_delta(rng);
    let current = xc.snap_up(rng).get_delta();

    for (current_iteration, ((i, t), _)) in thresh.into_iter().enumerate().cartesian_product(0..i_tilde as u64).enumerate() {
        let current_iteration = (current_iteration + 1) as f64;
        let total_iterations = ((i_tilde as u64) * (i_tilde as u64)) as f64;
        let k = points_grid.ordered_points.iter().map(|points| {
            let v = (points.len()-1)/2;
            v as f64 * (total_iterations-current_iteration)/total_iterations + current_iteration/total_iterations
        });
        let k: Vec<_> = k.map(|v| v as usize).collect();
        let mc: usize = (2.0 + current_iteration/total_iterations*(points_grid.d as f64 - 2.0)) as usize;
        let xc_plus = xc.generate_neighbor_delta(&k, mc, rng);
        let fxc = xc_plus.snap_up(rng).get_delta();
    }
}