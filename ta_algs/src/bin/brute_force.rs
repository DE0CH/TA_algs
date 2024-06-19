use ta_algs::entrance::get_raw_points;
use ta_algs::common::PointsGrid;
use ta_algs::common::PointIndex;
use ta_algs::common::Point;
use std::iter::repeat_with;
use itertools::Itertools;
use ordered_float::NotNan;

fn count_open(point: &Point) -> u64 {
    let points_grid = point.points_grid;
    (1..points_grid.n + 1).map(|i| {
        (0..points_grid.d).all(|id| {
            let x = points_grid.points[id][i];
            let y = points_grid.points[id][point.coord[id].0];
            x < y
        })
    }).map(|x| {
        let x: u64 = if x {1} else {0};
        return x
    }).sum()
}

fn count_closed(point: &Point) -> u64 {
    let points_grid = point.points_grid;
    (1..points_grid.n + 1).map(|i| {
        (0..points_grid.d).all(|id| {
            let x = points_grid.points[id][i];
            let y = points_grid.points[id][point.coord[id].0];
            x <= y
        })
    }).map(|x| {
        let x: u64 = if x {1} else {0};
        return x
    }).sum()
}

fn volume(point: &Point) -> NotNan<f64> {
    let points_grid = point.points_grid;
    (0..points_grid.d).map(|id| {
        let x = points_grid.points[id][point.coord[id].0];
        x
    }).product()

}


fn main() {
    let (points, _, _) = get_raw_points();
    let points_grid = PointsGrid::new(points);
    let grid_points = repeat_with(|| (1..points_grid.n + 2).map(|i| PointIndex(i))).take(points_grid.d).multi_cartesian_product();



    let (delta, bar_delta): (Vec<_>, Vec<_>) = grid_points.map(|point| {
        let point = Point::new(point, &points_grid);
        let delta = point.get_delta();
        let bardelta = point.get_bardelta();
        (delta, bardelta)
    }).unzip();

    println!("delta:    {}", delta.into_iter().max().unwrap());
    println!("bardelta: {}", bar_delta.into_iter().max().unwrap());


    let grid_points = repeat_with(|| (1..points_grid.n + 2).map(|i| PointIndex(i))).take(points_grid.d).multi_cartesian_product();

    let (delta, bar_delta): (Vec<_>, Vec<_>) = grid_points.map(|point| {
        let point = Point::new(point, &points_grid);
        let delta = *volume(&point) - (count_open(&point) as f64 / points_grid.n as f64);
        let bardelta = (count_closed(&point) as f64 / points_grid.n as f64) - *volume(&point);
        (NotNan::new(delta).unwrap(), NotNan::new(bardelta).unwrap())
    }).unzip();


    println!("delta:    {}", delta.into_iter().max().unwrap());
    println!("bardelta: {}", bar_delta.into_iter().max().unwrap());

}
