use std::clone::Clone;

use ordered_float::NotNan;
use fastrand::Rng;

const ERROR_MSG_TOO_LOW: &'static str = "Smaller than the smallest point in the grid. Since the grid should contains 0, it means either the point is smaller than 0 or the grid is not properly constructed";
const ERROR_MSG_TOO_HIGH: &'static str = "Larger than the largest point in the grid. Since the grid should contains 1, it means either the point is bigger than 1 or the grid is not properly constructed";

pub enum PM {
    Add,
    Subtract,
}

pub struct PointsGrid {
    pub d: usize,
    pub n: usize,
    pub raw_points: Vec<Vec<NotNan<f64>>>,
    pub coord: Vec<Vec<NotNan<f64>>>,
    pub point_index: Vec<Vec<usize>>
}

pub struct Point<'a> {
    coord: Vec<NotNan<f64>>,
    points_grid: &'a PointsGrid,
}

pub struct PointIndex<'a> {
    coord: Vec<usize>,
    points_grid: &'a PointsGrid,
}

impl<'a> Clone for Point<'a> {
    fn clone(&self) -> Self {
        Point {
            coord: self.coord.clone(),
            points_grid: self.points_grid,
        }
    }
}

impl<'a> Clone for PointIndex<'a> {
    fn clone(&self) -> Self {
        PointIndex {
            coord: self.coord.clone(),
            points_grid: self.points_grid,
        }
    }
}

pub struct K(pub Vec<usize>);

impl<'a> Point<'a> {
    pub fn new(points_grid: &'a PointsGrid) -> Self {
        Point {
            coord: Vec::with_capacity(points_grid.d),
            points_grid,
        }
    }

    pub fn new_random(points_grid: &'a PointsGrid, rng: &mut Rng) -> Self {
        Point {
            coord: (0..points_grid.d).map(|_| points_grid.wrap_with_pd_zo(rng.f64())).collect(),
            points_grid: &points_grid,
        }
    }

    pub fn round_point_up(&self) -> PointIndex<'a> {
        PointIndex {
            coord: self.coord.iter().enumerate().map(|(id, v)| self.get_index_up(id, *v).expect(ERROR_MSG_TOO_HIGH)).collect(),
            points_grid: self.points_grid,
        }
    }

    pub fn round_point_down(&self) -> PointIndex<'a> {
        PointIndex {
            coord: self.coord.iter().enumerate().map(|(id, v)| self.get_index_down(id, *v).expect(ERROR_MSG_TOO_LOW)).collect(),
            points_grid: self.points_grid,
        }
    }

    pub fn round_point_extradown(&self) -> PointIndex<'a> {
        //TODO: check if this is correct according to the paper
        PointIndex {
            coord: self.coord.iter().enumerate().map(|(id, v), | match self.get_index_down(id, *v).expect(ERROR_MSG_TOO_LOW) {
                0 => self.points_grid.coord[id].len() - 1,
                index => index,
            }).collect(),
            points_grid: &self.points_grid,
        }
    }

    pub fn get_index_up(&self, i: usize, v: NotNan<f64>) -> Option<usize> {
        let search_vector = &self.points_grid.coord[i];
        if v < *search_vector.get(0)? {
            return Some(0);
        }
        if v > *search_vector.last()? {
            return None
        }
        let mut left: usize = 0;
        let mut right: usize = search_vector.len();
        while left + 1 < right {
            let mid = (left + right) / 2;
            if search_vector[mid] < v {
                left = mid;
            } else {
                right = mid
            }
        }
        Some(left + 1)
    }

    pub fn get_index_down(&self, i: usize, v: NotNan<f64>) -> Option<usize> {
        let search_vector = &self.points_grid.coord[i];
        if v < *search_vector.get(0)? {
            return None;
        }
        let mut left: usize = 0;
        let mut right: usize = search_vector.len();
        while left + 1 < right {
            let mid = (left + right) / 2;
            if search_vector[mid] <= v {
                left = mid;
            } else {
                right = mid;
            }
        }
        Some(left)
    }

}

impl<'a> PointIndex<'a> {
    pub fn point(&self) -> Point<'a> {
        Point {
            coord: self.coord.iter().enumerate().map(|(id, j)| self.points_grid.coord[id][*j]).collect(),
            points_grid: self.points_grid,
        }
    }

    pub fn generate_neighbor(&self, k: K, mc: usize, rng: &mut fastrand::Rng) -> (PointIndex, PointIndex, PointIndex){
        let d = self.points_grid.d;
        let mut point = self.point();
        let rand_coord = rng.choose_multiple(0..d, mc);
        
        for id in rand_coord.into_iter() {
            let diff = k.0[id];
            let lower_bound = self.get_clipped(id, diff, PM::Subtract);
            let upper_bound = self.get_clipped(id, diff, PM::Add);
            let s = rng.f64();
            point.coord[id] = self.points_grid.wrap_with_pd(lower_bound, upper_bound, s);
        }

        (point.round_point_up(), point.round_point_down(), point.round_point_extradown())
    }

    pub fn generate_neighbor_delta(&self, k: K, mc: usize, rng: &mut fastrand::Rng) -> PointIndex {
        let d = self.points_grid.d;
        let mut point = self.point();
        let rand_coord = rng.choose_multiple(0..d, mc);
        
        for id in rand_coord.into_iter() {
            let diff = k.0[id];
            let lower_bound = self.get_clipped(id, diff, PM::Subtract);
            let upper_bound = self.get_clipped(id, diff, PM::Add);
            let s = rng.f64();
            point.coord[id] = self.points_grid.wrap_with_pd(lower_bound, upper_bound, s);
        }

        point.round_point_up()
    }

    pub fn get_clipped(&self, id: usize, diff: usize, pm: PM) -> NotNan<f64> {
        self.points_grid.get_clipped(id, self.coord[id], diff, pm)
    }
    
}

impl<'a> PointsGrid {
    
    pub fn new(raw_points: Vec<Vec<NotNan<f64>>>) -> Self {
        let n = raw_points.len();
        let d = raw_points[0].len();
        let n_dimensions = d;
        let n_points = n;
    
        // let coordinate: Vec<usize> = (0..d).collect();
        
        let (coord, point_index): (Vec<Vec<NotNan<f64>>>, Vec<Vec<usize>>) = (0..n_dimensions).map(|i| -> (Vec<NotNan<f64>>, Vec<usize>) {
    
            let zero = NotNan::new(0.0f64).unwrap();
            let zero = (0usize, zero);
            let one = NotNan::new(1.0f64).unwrap();
            let one = (n_points+1, one);
            let mut temp_coord: Vec<(usize, NotNan<f64>)> = std::iter::once(zero)
                .chain((0..n_points).map(|j| (j+1, raw_points[j][i])))
                .chain(std::iter::once(one))
                .collect();
    
            temp_coord.sort_by(|a, b| a.1.cmp(&b.1));
    
            let (index_map, temp_coord) = dedup_and_get_coord(temp_coord);
            
            let ans = temp_coord.into_iter().map(|a| a.1).collect();
            
            (ans, index_map)
    
        }).unzip();
    
        PointsGrid {
            n, 
            d,
            raw_points,
            coord,
            point_index,
        }
    }

    pub fn generate_xc(&'a self, rng: &mut Rng) -> (PointIndex<'a>, PointIndex<'a>, PointIndex<'a>) {
        let point = Point::new_random(&self, rng);
        (point.round_point_up(), point.round_point_down(), point.round_point_extradown())
    }

    pub fn generate_xc_delta(&'a self, rng: &mut Rng) -> PointIndex<'a> {
        let point = Point::new_random(&self, rng);
        point.round_point_up()
    }

    pub fn generate_xc_bardelta(&'a self, rng: &mut Rng) -> PointIndex<'a> {
        let point = Point::new_random(&self, rng);
        point.round_point_up()
    }

    pub fn get_clipped(&self, id: usize, i: usize, diff: usize, pm: PM) -> NotNan<f64> {
        //pm: true is plus, pm false is minus
        match pm {
            PM::Add => {
                match i.checked_add(diff) {
                    Some(x) => match self.coord[id].get(x) {
                        Some(v) => *v,
                        None => NotNan::new(1.0).unwrap(),
                    },
                    None => NotNan::new(1.0).unwrap(),
                }
            }
            PM::Subtract => {
                match i.checked_sub(diff) {
                    Some(x) => match self.coord[id].get(x) {
                        Some(v) => *v,
                        None => NotNan::new(1.0).unwrap(),
                    }
                    None => NotNan::new(1.0).unwrap(),
                }
            }
        }
    }

    pub fn wrap_with_pd_zo(&self, s: f64) -> NotNan<f64> {
        self.wrap_with_pd(NotNan::new(0.0).unwrap(), NotNan::new(1.0).unwrap(), s)
    }

    pub fn wrap_with_pd(&self, lower_bound: NotNan<f64>, upper_bound: NotNan<f64>, s: f64) -> NotNan<f64> {
        let d = self.d;
        let temp = (upper_bound.powi(d as i32) - lower_bound.powi(d as i32))*s + lower_bound.powi(d as i32);
        NotNan::new(temp.powf(1.0 / d as f64)).unwrap()
    }
}

// remove duplicate from indexed elements and return the permutation
pub fn dedup_and_get_coord<T: PartialEq>(mut stuff: Vec<(usize, T)>) -> (Vec<usize>, Vec<(usize, T)>) {
    let mut resulting_map: Vec<usize> = Vec::with_capacity(stuff.len());
    // TODO: check this actual does the correct thing
    stuff.dedup_by(|a, b| {
        match a.1 == b.1 {
            true => {
                resulting_map[b.0] = a.0;
                true
            },
            false => {
                resulting_map[b.0] = b.0;
                false
            }
        }
    });
    (resulting_map, stuff)
}

