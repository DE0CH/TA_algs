use std::clone::Clone;

use ordered_float::NotNan;
use fastrand::Rng;
use std::iter::{zip, repeat, repeat_with};
use itertools::izip;
use std::cmp::Ordering;
use std::cmp::Ord;
use crate::utils::into_iterator_with_permutation_unchecked;
use crate::utils::print_iterator;

const ERROR_MSG_TOO_LOW: &'static str = "Smaller than the smallest point in the grid. Since the grid should contains 0, it means either the point is smaller than 0 or the grid is not properly constructed";
const ERROR_MSG_TOO_HIGH: &'static str = "Larger than the largest point in the grid. Since the grid should contains 1, it means either the point is bigger than 1 or the grid is not properly constructed";

pub enum PM {
    Add,
    Subtract,
}

pub struct PointsGrid {
    pub d: usize,
    pub n: usize,
    pub points: Vec<Vec<NotNan<f64>>>,
    pub zeros: Vec<NotNan<f64>>,
    pub ones: Vec<NotNan<f64>>,
    pub ordered_points: Vec<Vec<NotNan<f64>>>,
    pub permutation: Vec<Vec<usize>>,
}

#[derive(Clone, Copy, Debug)]
pub enum Index {
    Zero,
    One,
    N(usize),
}

impl Ord for Index {
    fn cmp(&self, other: &Self) -> Ordering {
        match (self, other) {
            (Index::Zero, Index::Zero) => Ordering::Equal,
            (Index::Zero, Index::One) => Ordering::Less,
            (Index::Zero, Index::N(_)) => Ordering::Less,
            (Index::One, Index::Zero) => Ordering::Greater,
            (Index::One, Index::One) => Ordering::Equal,
            (Index::One, Index::N(_)) => Ordering::Less,
            (Index::N(_), Index::Zero) => Ordering::Greater,
            (Index::N(_), Index::One) => Ordering::Greater,
            (Index::N(x), Index::N(y)) => x.cmp(y),
        }
    
    }
}

impl PartialOrd for Index {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}

impl PartialEq for Index {
    fn eq(&self, other: &Self) -> bool {
        match (self, other) {
            (Index::Zero, Index::Zero) => true,
            (Index::One, Index::One) => true,
            (Index::N(x), Index::N(y)) => x == y,
            _ => false,
        }
    }
}

impl Eq for Index {}

impl Index {
    fn shift_index(&self, max: &usize) -> Option<usize> {
        match self {
            Index::Zero => Some(0),
            Index::N(x) => Some(x.checked_add(1))?,
            Index::One => Some(max.checked_add(1))?,
        }
    }
    fn add(&self, diff: &usize, max: &usize) -> Option<Index> {
        let shifted_index = self.shift_index(max)?;
        match shifted_index.checked_add(*diff)? {
            max if max == max + 1 => Some(Index::One),
            0 => Some(Index::Zero),
            x => Some(Index::N(x-1)),
        }
    }
    fn subtract(&self, diff: &usize, max: &usize) -> Option<Index> {
        let shifted_index = self.shift_index(max)?;
        match shifted_index.checked_sub(*diff)? {
            0 => Some(Index::Zero),
            max if max == max + 1 => Some(Index::One),
            x => Some(Index::N(x-1)),
        }

    }
}

#[derive(Clone, Copy)]
pub enum RoughIndex {
    Index(Index),
    Float(NotNan<f64>),
}

pub struct RoughPoint<'a> {
    pub coord: Vec<RoughIndex>,
    pub points_grid: &'a PointsGrid,
}

pub struct Point<'a> {
    pub coord: Vec<Index>,
    pub points_grid: &'a PointsGrid,
}

impl<'a> RoughPoint<'a> {
    pub fn new_random(points_grid: &'a PointsGrid, rng: & mut Rng) -> Self {
        RoughPoint {
            coord: repeat_with(|| RoughIndex::Float(NotNan::new(rng.f64()).unwrap())).take(points_grid.d).collect(),
            points_grid: &points_grid,
        }
    }
    pub fn get_index_up(&'a self) -> impl Iterator<Item = Option<Index>> + 'a {
        self.points_grid.get_index_up(self.coord.iter().map(|x| *x))
    }

    pub fn get_index_down(&'a self) -> impl Iterator<Item = Option<Index>> + 'a {
        self.points_grid.get_index_down(self.coord.iter().map(|x| *x))
    }

    fn round_point<'b>(rounded_index: impl Iterator<Item = Option<Index>> + 'b, non_choice: impl Iterator<Item = Index> + 'b) -> impl Iterator<Item = Index> + 'b
    where
        'a: 'b
    {
        izip!(rounded_index, non_choice).map(|(rounded_index, non_choice)| {
            match rounded_index {
                Some(x) => x,
                None => non_choice,
            }
        })
    }

    pub fn round_point_up(&self) -> Point<'a> {
        let points_grid = self.points_grid;
        Point {
            coord: Self::round_point(self.get_index_up(), repeat(Index::One)).collect(),
            points_grid,
        }
    }

    pub fn round_point_down(&self) -> Point<'a> {
        let points_grid = self.points_grid;
        Point {
            coord: Self::round_point(self.get_index_up(), repeat(Index::Zero)).collect(),
            points_grid
        }
    }

    pub fn round_point_extradown(&self) -> Point<'a> {
        let points_grid = self.points_grid;
        Point {
            coord: Self::round_point(self.get_index_up(), repeat(Index::N(points_grid.n))).collect(),
            points_grid
        }
    }
}

impl<'a> Point<'a> {

    pub fn generate_neighbor_point<'b>(&self, k: &'b [usize], mc: usize, rng: &'b mut Rng) -> RoughPoint<'a>
    where
        'a: 'b
    {
        let d = self.points_grid.d;
        let rand_coord = rng.choose_multiple(0..d, mc);
        let rand_coord = subset_to_mask(&rand_coord, d);
        let lower_bound = self.get_clipped(k.iter(), PM::Subtract);
        let upper_bound = self.get_clipped(k.iter(), PM::Add);
        let s = repeat_with(|| rng.f64()).take(d).collect::<Vec<f64>>();
        
        let coord = izip!(rand_coord, s, lower_bound, upper_bound, self.coord.iter()).map(|(rand_coord, s, lower_bound, upper_bound, i)| {
            if rand_coord {
                RoughIndex::Float(self.points_grid.warp_with_pd(*lower_bound, *upper_bound, s))
            } else {
                RoughIndex::Index(*i)
            }
        });
        RoughPoint {
            coord: coord.collect(),
            points_grid: self.points_grid,
        }
    }

    pub fn get_delta(&self) -> NotNan<f64> {
        let n = self.points_grid.n;
        let op = self.count_open();
        let vol = self.volume();
        vol - NotNan::new(op as f64 / n as f64).unwrap()
    }

    pub fn get_bardelta(&self) -> NotNan<f64> {
        let n = self.points_grid.n;
        let cl = self.count_closed();
        let vol = self.volume();
        NotNan::new(cl as f64 / n as f64).unwrap() - vol
    }

    pub fn volume(&self) -> NotNan<f64> {
        izip!(self.points_grid.zeros.iter(), self.points_grid.ones.iter(), self.points_grid.points.iter(), self.coord.iter()).map(|(zero, one, point, v)| match v {
            Index::Zero => zero,
            Index::One => one,
            Index::N(i) => &point[*i],
        }).map(|x| *x).reduce(|a, b| a * b).unwrap()
    }

    // is x in open box bounded by self?, since we are counting only points, we don't use Index for x, but just usize
    pub fn open(&self, x: impl Iterator<Item = &'a NotNan<f64>>) -> bool {
        izip!(self.points_grid.zeros.iter(), self.points_grid.ones.iter(), x, self.coord.iter(), self.points_grid.points.iter()).all(|(zero, one, x, c, points)| {
            let c: &NotNan<f64> = match c {
                Index::Zero => zero,
                Index::One => one,
                Index::N(i) => &points[*i],
            };
            *x < *c
        })
    }

    pub fn closed(&self, x: impl Iterator<Item = &'a NotNan<f64>>) -> bool {
        izip!(self.points_grid.zeros.iter(), self.points_grid.ones.iter(), x, self.coord.iter(), self.points_grid.points.iter()).all(|(zero, one, x, c, points)| {
            let c: &NotNan<f64> = match c {
                Index::Zero => zero,
                Index::One => one,
                Index::N(i) => &points[*i],
            };
            *x <= *c
        })
    }

    pub fn count_open(&self) -> usize {
        let xs = self.points_grid.point_set();
        xs.map(|(_, x)| {
            if self.open(x) {1} else {0} 
        }).fold(0, |a, b| a + b)
    }

    pub fn count_closed(&self) -> usize {
        let xs = self.points_grid.point_set();
        xs.map(|(_, x)| {
            if self.closed(x) {1} else {0}
        }).fold(0, |a, b| a + b)
    }

    pub fn get_clipped<'b>(&'a self, diff: impl Iterator<Item = &'b usize> + 'b, pm: PM) -> impl Iterator<Item = &NotNan<f64>> + 'b 
    where 
        'a: 'b
    {
        self.points_grid.get_clipped(self.coord.iter(), diff, pm)
    }
    
    pub fn snap_up(&self, rng: &mut Rng) -> Point<'a> {
        let order = self.points_grid.get_permutation(rng);
        let reverse_order = get_reveres_permutation(&order);
        let new_box = self.points_grid.new_max_index();
        let xs = self.points_grid.point_set();
        let xs2 = self.points_grid.point_set();
        let xs3 = self.points_grid.point_set_with_permutation(&order);
        izip!(xs, xs2, xs3).fold(new_box, |new_box, ((_, x1), (_, x2), (x3i, x3))| {
            if !self.open(x1) && new_box.open(x2) {
                Point {
                    coord: reverse_iterator_with_permutation({
                        let yp = unsafe {into_iterator_with_permutation_unchecked(self.to_float().collect::<Vec<_>>(), &order)};
                        let yp = yp.collect::<Vec<_>>();
                        let yp = yp.into_iter();
                        let yp_sn = iterator_with_permutation(&new_box.coord, &order);
                        izip!(repeat(x3i), x3, yp, yp_sn)
                        .scan(false, |found, (xi, x, yp, yp_sn)| {
                            if !*found && *x >= *yp {
                                *found = true;
                                Some(Index::N(xi))
                            } else {
                                Some(*yp_sn)
                            } 
                        })
                    }, &reverse_order),
                    points_grid: self.points_grid,
                }
            } else {
                new_box
            }
        })
    }

    pub fn snap_down(&self) -> Point<'a> {
        let xs = self.points_grid.point_set();
        let xs2 = self.points_grid.point_set();
        let point = self.points_grid.new_origin_index();

        zip(xs, xs2).fold(point, |point, ((_, x1), (x2i, x2))| {
                if self.closed(x1) {
                    Point { 
                        coord: izip!(point.coord.iter(), point.to_float(), repeat(x2i), x2).map(|(pi, p, xi, x)| {
                            if *p > *x {
                                Index::N(xi)
                            } else {
                                *pi
                            }
                        }).collect(),
                        points_grid: self.points_grid,
                    }
                } else {
                    point
                }
        })
    }

    fn to_float(&'a self) -> impl Iterator<Item = &'a NotNan<f64>> + 'a {
        izip!(self.points_grid.zeros.iter(), self.points_grid.ones.iter(), self.points_grid.points.iter(), self.coord.iter()).map(|(zero, one, points, i)| {
            match i {
                Index::Zero => zero,
                Index::One => one,
                Index::N(i) => &points[*i],
            }
        })
    }
}

impl<'a> PointsGrid {
    
    pub fn new(points: Vec<Vec<NotNan<f64>>>) -> Self {
        let d = points.len();
        let n = points[0].len();
        let (ordered_points, permutation) = points.iter().map(|points| {
            let mut indexed_points: Vec<_> = points.into_iter().map(|x| *x).enumerate().collect();
            indexed_points.sort_by_key(|x| x.1);
            let (permutation, ordered_points) = dedup_and_get_coord(indexed_points);
            (ordered_points, permutation)
        }).unzip();
        PointsGrid {
            n, 
            d,
            points,
            ordered_points,
            permutation,
            zeros: vec![NotNan::new(0.0).unwrap(); d],
            ones: vec![NotNan::new(1.0).unwrap(); d],
        }
    }

    pub fn new_origin_index(&'a self) -> Point<'a> {
        Point {
            coord: repeat(Index::Zero).take(self.d).collect(),
            points_grid: self,
        }
    }

    pub fn new_max_index(&'a self) -> Point<'a> {
        Point {
            coord: repeat(Index::One).take(self.d).collect(),
            points_grid: self,
        }
    }

    pub fn get_permutation(&self, rng: &mut Rng) -> Vec<usize> {
        let mut order: Vec<_> = (0..self.d).collect();
        rng.shuffle(&mut order);
        order
    }

    pub fn generate_random_point(&'a self, rng: &mut Rng) -> RoughPoint<'a> {
        RoughPoint::new_random(&self, rng)
    }

    pub fn get_clipped<'b>(&'a self, i: impl Iterator<Item = &'a Index> + 'b, diff: impl Iterator<Item = &'b usize> + 'b, pm: PM) -> impl Iterator<Item = &'a NotNan<f64>> + 'b 
    where 
        'a: 'b 
    {
        izip!(self.zeros.iter(), self.ones.iter(), i, diff, self.points.iter()).map(move |(zero, one, i, diff, coord)| {
            match pm {
                PM::Add => {
                    match i.add(diff, &coord.len()) {
                        Some(x) => match x {
                            Index::One => one,
                            Index::Zero => zero,
                            Index::N(i) => &coord[i],
                        },
                        None => one,
                    }
                },
                PM::Subtract => {
                    match i.subtract(diff, &coord.len()) {
                        Some(x) => match x {
                            Index::One => one, 
                            Index::Zero => zero,
                            Index::N(i) => &coord[i],
                        }
                        None => zero
                    }
                }
            }
        })
    }

    pub fn point_set(&self) -> impl Iterator<Item = (usize, impl Iterator<Item = &NotNan<f64>>)> {
        (0..self.n).map(move |i| (i, (0..self.d).map(move |id| {
            &self.points[id][i]
        })))
    }

    pub fn point_set_with_permutation<'b>(&'a self, permutation: &'b [usize]) -> impl Iterator<Item = (usize, impl Iterator<Item = &'a NotNan<f64>> + 'b)> + 'b where 'a: 'b{
        let ans = (0..self.n).map(move |i| {
            (i, (0..self.d).map(move |id| {
                &self.points[permutation[id]][i]
            }))
        });
        ans
    }

    pub fn warp_with_pd_zo(&self, s: f64) -> NotNan<f64> {
        self.warp_with_pd(NotNan::new(0.0).unwrap(), NotNan::new(1.0).unwrap(), s)
    }

    pub fn warp_with_pd(&self, lower_bound: NotNan<f64>, upper_bound: NotNan<f64>, s: f64) -> NotNan<f64> {
        let d = self.d;
        let temp = (upper_bound.powi(d as i32) - lower_bound.powi(d as i32))*s + lower_bound.powi(d as i32);
        NotNan::new(temp.powf(1.0 / d as f64)).unwrap()
    }
  
    pub fn get_index_up<'b>(&'a self, v: impl Iterator<Item = RoughIndex> + 'b) -> impl Iterator<Item = Option<Index>> + 'b
    where
        'a: 'b
    {
        izip!(self.points.iter(), self.permutation.iter(), v).map(|(search_vector, permutation, v)| {
            match v {
                RoughIndex::Float(v) => {
                    match get_index_up(search_vector, &v) {
                        Some(x) => Some(Index::N(permutation[x])),
                        None => None,
                    }
                },
                RoughIndex::Index(x) => Some(x)
            }
        })
    }

    pub fn get_index_down<'b>(&'a self, v: impl Iterator<Item = RoughIndex> + 'b) -> impl Iterator<Item = Option<Index>> + 'b
    where
        'a: 'b
    {
        izip!(self.points.iter(), self.permutation.iter(), v).map(|(search_vector, permutation, v)| {
            match v {
                RoughIndex::Float(v) => {
                    match get_index_down(search_vector, &v) {
                        Some(x) => Some(Index::N(permutation[x])),
                        None => None,
                    }
                },
                RoughIndex::Index(x) => Some(x)
            }
        })
    }

}

// remove duplicate from indexed elements and return the permutation
pub fn dedup_and_get_coord<T: PartialEq>(stuff: Vec<(usize, T)>) -> (Vec<usize>, Vec<T>) {
    let mut order = Vec::<usize>::with_capacity(stuff.len());
    let mut just_things = Vec::<T>::with_capacity(stuff.len());
    stuff.into_iter().for_each(|(ii, vv)| {
        if let Some(v) = just_things.last() {
            if vv != *v {
                order.push(ii);
                just_things.push(vv);
            }
        } else {
            order.push(ii);
            just_things.push(vv);
        }
    });
    (order, just_things)
}

fn iterator_with_permutation<'b, T>(stuff: &'b [T], permutation: &'b [usize]) -> impl Iterator<Item = &'b T> + 'b {
    let iterator = permutation.iter().map(|i| &stuff[*i]);
    iterator
}

fn get_reveres_permutation(permutation: &[usize]) -> Vec<usize> {
    let mut reverse_permutation = Vec::with_capacity(permutation.len());
    unsafe {reverse_permutation.set_len(permutation.len())}
    for (i, j) in permutation.iter().enumerate() {
        reverse_permutation[*j] = i;
    }
    reverse_permutation
}

fn reverse_iterator_with_permutation<T>(iterator: impl Iterator<Item = T>, reverse_permutation: &[usize]) -> Vec<T> {
    let mut ans = Vec::with_capacity(reverse_permutation.len());
    unsafe {ans.set_len(reverse_permutation.len())}
    for (i, v) in iterator.enumerate() {
        ans[reverse_permutation[i]] = v;
    }
    ans
}

// TODO: test the algorithms 
fn get_index_up_with_cmp<T, P>(search_vector: &[T], mut cmp: P) -> Option<usize>
where 
    P: FnMut(&T) -> Ordering
{
    if cmp(search_vector.get(0)?) == Ordering::Greater {
        return Some(0);
    }
    if cmp(search_vector.last()?) == Ordering::Less {
        return None;
    }
    let mut left: usize = 0;
    let mut right: usize = search_vector.len();
    while left + 1 < right {
        let mid = (left + right) / 2;
        if cmp(&search_vector[mid]) == Ordering::Less {
            left = mid;
        } else {
            right = mid;
        }
    }
    // Because of the guard at the beginning, it's guaranteed that left + 1 is within bounds. 
    Some(left + 1)
}

fn get_index_up<T: Ord>(search_vector: &[T], v: &T) -> Option<usize> {
    get_index_up_with_cmp(search_vector, |other| other.cmp(v))
}

fn get_index_down_with_cmp<T, P>(search_vector: &[T], mut cmp: P) -> Option<usize>
where
    P: FnMut(&T) -> Ordering
{
    if cmp(search_vector.get(0)?) == Ordering::Greater {
        return None;
    }
    let mut left: usize = 0;
    let mut right: usize = search_vector.len();
    while left + 1 < right {
        let mid = (left + right) / 2;
        match cmp(&search_vector[mid]) {
            Ordering::Less | Ordering::Equal => left = mid,
            Ordering::Greater => right = mid,
        };
    }

    Some(left)
}

fn get_index_down<T: Ord>(search_vector: &[T], v: &T) -> Option<usize> {
    get_index_down_with_cmp(search_vector, |other| other.cmp(v))
}

fn subset_to_mask(subset: &[usize], n: usize) -> Vec<bool> {
    let mut mask = vec![false; n];
    for i in subset {
        mask[*i] = true;
    }
    mask
}