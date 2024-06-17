use ordered_float::NotNan;
use ta_algs::common::{Index, PointsGrid};
use fastrand::Rng;
use itertools::Itertools;
use core::cmp::max;
use ta_algs::entrance;

pub fn trial(points_grid: PointsGrid, iterations: u64, rng: &mut Rng) -> (Vec<Index>, NotNan<f64>, NotNan<f64>){
    let i_tilde = (iterations as f64).sqrt() as usize;
    let mut thresh: Vec<_> = (0..i_tilde as usize).map(|i| {
        let current_iteration = (i+1) as f64;
        let i_tilde = i_tilde as f64;
        let k = points_grid.ordered_points.iter().map(|points| {
            ((points.len()-1) as f64)/2.0 * (i_tilde-current_iteration)/i_tilde + current_iteration/i_tilde
        });
        let k: Vec<_> = k.map(|v| v as usize).collect();
        let mc: usize = (2.0 + current_iteration/i_tilde*(points_grid.d as f64 - 2.0)) as usize;
        let x = points_grid.generate_random_point(rng);
        let xm = x.round_point_down();
        let xm_sn = xm.snap_down();
        let y = xm.generate_neighbor_point(&k, mc, rng);
        let ym = y.round_point_down();
        let ym_sn = ym.snap_down();
        let ymm = y.round_point_extradown();
        let ymm_sn = ymm.snap_down();

        NotNan::new(-(xm_sn.get_bardelta() - max(ym_sn.get_bardelta(), ymm_sn.get_bardelta())).abs()).unwrap()
    }).collect();
    thresh.sort();

    let x = points_grid.generate_random_point(rng);
    let xc = x.round_point_up();
    let global = xc.snap_up(rng).get_bardelta();
    let current = xc.snap_up(rng).get_bardelta();
    let big_i = (i_tilde * i_tilde) as f64;
    let ans = thresh.into_iter().cartesian_product(0..i_tilde).enumerate().fold((xc, global, current), |(xc, global, current), (current_iteration, (tt, _))|
        {
            let current_iteration = (current_iteration + 1) as f64;
            let k = points_grid.ordered_points.iter().map(|points| {
                ((points.len()-1) as f64)/(2.0 as f64) * (big_i-current_iteration)/big_i + current_iteration/big_i
            });
            let k: Vec<_> = k.map(|v| v as usize).collect();
            let mc: usize = (2.0 + current_iteration/big_i*(points_grid.d as f64 - 2.0)) as usize;
            let y = xc.generate_neighbor_point(&k, mc, rng);
            let ym = y.round_point_down();
            let ym_sn = ym.snap_down();
            let ymm = y.round_point_extradown();
            let ymm_sn = ymm.snap_down();
            let delta = max(ym_sn.get_bardelta(), ymm_sn.get_bardelta());
            let ys = if ym_sn.get_bardelta() >= ymm_sn.get_bardelta() {ym} else {ymm};
            let global = max(global, delta);
            let change = delta - current;
            if change >= tt {
                (ys, global, delta)
            } else {
                (xc, global, current)
            }
        }
    );
    let (xc, global, current) = ans;
    let xc = xc.coord;
    (xc, global, current)
}

fn main() {
    let (points, seed, iterations) = entrance::get_raw_points();
    let mut rng = fastrand::Rng::with_seed(seed);
    let points2 = points.clone();
    let points = PointsGrid::new(points);
    let (a, b, c) = trial(points, iterations, &mut rng);
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
