use std::clone::Clone;
use std::iter::once;

use ordered_float::NotNan;
use fastrand::Rng;
use std::iter::{zip, repeat, repeat_with};
use itertools::izip;
use std::cmp::Ordering;
use std::cmp::Ord;
use itertools::Itertools;

const ERROR_MSG_TOO_LOW: &'static str = "Smaller than the smallest point in the grid. Since the grid should contains 0, it means either the point is smaller than 0 or the grid is not properly constructed";
const ERROR_MSG_TOO_HIGH: &'static str = "Larger than the largest point in the grid. Since the grid should contains 1, it means either the point is bigger than 1 or the grid is not properly constructed";

pub enum PM {
    Add,
    Subtract,
}

#[derive(Clone, Copy, Debug, Eq, PartialEq, Ord, PartialOrd)]
pub struct PointIndex(pub usize);
#[derive(Clone, Copy, Debug, Eq, PartialEq, Ord, PartialOrd)]
pub struct PointOrder(pub usize);

pub struct PointsGrid {
    pub d: usize,
    pub n: usize,
    pub points: Vec<Vec<NotNan<f64>>>,
    pub ordered_points: Vec<Vec<NotNan<f64>>>,
    pub reverse_order: Vec<Vec<PointIndex>>,
    pub order: Vec<Vec<PointOrder>>,
}

#[derive(Clone, Copy)]
pub enum RoughIndex {
    Index(PointIndex),
    Float(NotNan<f64>),
}

pub struct RoughPoint<'a> {
    pub coord: Vec<RoughIndex>,
    pub points_grid: &'a PointsGrid,
}

pub struct Point<'a> {
    pub coord: Vec<PointIndex>,
    pub points_grid: &'a PointsGrid,
}

impl<'a> RoughPoint<'a> {
    pub fn new_random(points_grid: &'a PointsGrid, rng: & mut Rng) -> Self {
        RoughPoint {
            coord: repeat_with(|| {
                let x = NotNan::new(rng.f64()).unwrap();
                let x = points_grid.warp_with_pd_zo(x);
                RoughIndex::Float(x)
            }).take(points_grid.d).collect(),
            points_grid: &points_grid,
        }
    }

    pub fn round_up(&self) -> Point<'a> {
        let coord = zip(self.coord.iter(), self.points_grid.ordered_points.iter()).map(|(x, ordered_points)| {
            match x {
                RoughIndex::Index(i) => *i,
                RoughIndex::Float(f) => {
                    let ans = get_index_up(&ordered_points, &f).expect(ERROR_MSG_TOO_HIGH);
                    PointIndex(ans)
                }
            }
        }).collect();
        Point::new(
            coord,
            self.points_grid,
        )
    }

    pub fn round_down(&self) -> Point<'a> {
        let coord = zip(self.coord.iter(), self.points_grid.ordered_points.iter()).map(|(x, ordered_points)| {
            match x {
                RoughIndex::Index(i) => *i,
                RoughIndex::Float(f) => {
                    let ans = get_index_down(&ordered_points, &f).expect(ERROR_MSG_TOO_LOW);
                    PointIndex(ans)
                }
            }
        }).collect();
        Point::new(
            coord,
            self.points_grid,
        )
    }
}


impl<'a> Point<'a> {

    pub fn new(coord: Vec<PointIndex>, points_grid: &'a PointsGrid) -> Self {
        Point {
            coord,
            points_grid,
        }
    }

    pub fn new_from_order(coord: impl Iterator<Item = PointOrder>, points_grid: &'a PointsGrid) -> Self {
        let coord = zip(coord, &points_grid.reverse_order).map(|(x, reverse_order)| reverse_order[x.0]).collect();
        Point {
            coord,
            points_grid,
        }
    }

    pub fn generate_neighbor_point<'b>(&self, k: &'b [usize], mc: usize, rng: &'b mut Rng) -> RoughPoint<'a>
    where
        'a: 'b
    {
        let d = self.points_grid.d;
        let rand_coord = rng.choose_multiple(0..d, mc);
        let rand_coord = subset_to_mask(&rand_coord, d);
        let lower_bound = self.get_clipped(k.iter(), PM::Subtract);
        let upper_bound = self.get_clipped(k.iter(), PM::Add);
        let s = repeat_with(|| NotNan::new(rng.f64()).unwrap()).take(d).collect::<Vec<NotNan<f64>>>();

        let coord = izip!(rand_coord, s, lower_bound, upper_bound, self.coord.iter(), self.points_grid.points.iter()).map(|(rand_coord, s, lower_bound, upper_bound, i, points)| {
            if rand_coord {
                let lower_bound = points[lower_bound.0];
                let upper_bound = points[upper_bound.0];
                RoughIndex::Float(self.points_grid.warp_with_pd(lower_bound, upper_bound, s))
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
        let ans = NotNan::new(cl as f64 / n as f64).unwrap() - vol;
        ans
    }

    pub fn float(&'a self) -> impl Iterator<Item = NotNan<f64>> + 'a {
        let points = self.points_grid.points.iter();
        let coord = self.coord.iter();
        zip(points, coord).map(|(points, i)| points[i.0])
    }

    pub fn volume(&self) -> NotNan<f64> {
        let p = self.float();
        let one = NotNan::new(1.0).unwrap();
        p.fold(one, |a, b| a * b)
    }

    pub fn get_order(&'a self) -> impl Iterator<Item = PointOrder> + 'a {
        let order = self.points_grid.order.iter();
        let coord = self.coord.iter();
        zip(order, coord).map(|(order, i)| order[i.0])
    }

    // is x in open box bounded by self, since we are counting only points, we don't use Index for x, but just usize
    pub fn open(&self, x: PointIndex) -> bool {
        let self_ = self.get_order();
        let x = self.points_grid.point_order(x);
        zip(self_, x).all(|(s, x)| {
            x.0 < s.0
        })
    }

    pub fn closed(&self, x: PointIndex) -> bool {
        let self_ = self.get_order();
        let x = self.points_grid.point_order(x);
        zip(self_, x).all(|(s, x)| {
            x.0 <= s.0
        })
    }

    pub fn count_open(&self) -> usize {
        let xs = self.points_grid.point_set();
        xs.map(|x| {
            if self.open(x) {1} else {0}
        }).fold(0, |a, b| a + b)
    }

    pub fn count_closed(&self) -> usize {
        let xs = self.points_grid.point_set();
        xs.map(|x| {
            if self.closed(x) {1} else {0}
        }).fold(0, |a, b| a + b)
    }

    pub fn get_clipped<'b>(&'a self, diff: impl Iterator<Item = &'b usize> + 'b, pm: PM) -> impl Iterator<Item = PointIndex> + 'b
    where
        'a: 'b
    {
        let zero = self.points_grid.zero();
        let one = self.points_grid.one();
        let self_ = self.get_order();
        let ans = izip!(self_, diff, &self.points_grid.reverse_order).map(move |(x, diff, reverse_order)| {
            let x = x.0;
            let diff = *diff;
            let ans = match pm {
                PM::Add => x.checked_add(diff),
                PM::Subtract => x.checked_sub(diff),
            };
            let ans = match ans {
                Some(ans) => {
                    if ans == 0 {zero} else if ans >= reverse_order.len() {one} else {reverse_order[ans]}
                }
                None => {match pm{ PM::Add => one, PM::Subtract => zero}}
            };
            ans
        });
        ans
    }

    pub fn snap_up(&self, rng: &mut Rng) -> Point<'a> {
        let coordinate_order = self.points_grid.get_permutation(rng);
        let reverse_order: Vec<_> = get_reveres_permutation(&coordinate_order).into_iter().map(|x| x.unwrap()).collect();
        let xs1 = self.points_grid.point_set();
        let xs2 = self.points_grid.point_set();
        let xs3 = self.points_grid.point_set();

        let mut new_box = self.points_grid.new_max_index();

        for (x1, x2, x3) in izip!(xs1, xs2, xs3) {
            if !self.open(x1) && new_box.open(x2) {
                let yp_point_order: Vec<_> = new_box.get_order().collect();
                let yp_sn = iterator_with_permutation(&yp_point_order, &coordinate_order);
                let yp = iterator_with_permutation(&yp_point_order, &coordinate_order);
                let x3_point_order: Vec<_> = self.points_grid.point_order(x3).collect();
                let x = iterator_with_permutation(&x3_point_order, &coordinate_order);

                // OPTIMIZATION: update in place, though I suspect it will already be optimized by the compiler
                let new = izip!(x, yp, yp_sn).scan(false, |found, (x, yp, yp_sn)| {
                    if !*found && x >= yp {
                        *found = true;
                        Some(*x)
                    } else {
                        Some(*yp_sn)
                    }
                });

                // OPTIMIZATION: don't do this every time
                let new = reverse_iterator_with_permutation(new, &reverse_order).into_iter().map(|x| x.unwrap());
                new_box = Point::new_from_order(new, self.points_grid);
            }
        }

        new_box
    }

    pub fn snap_down(&self) -> Point<'a> {
        let xs1 = self.points_grid.point_set();
        let xs2 = self.points_grid.point_set();
        let mut new_box = self.points_grid.new_origin_index();

        for (x1, x2) in izip!(xs1, xs2) {
            if self.closed(x1) && !new_box.open(x2) {
                let x = self.points_grid.point_order(x1);
                let ans = izip!(new_box.get_order(), x).map(|(ym_sn, x)| {
                    if x > ym_sn {
                        x
                    } else {
                        ym_sn
                    }
                });
                new_box = Point::new_from_order(ans, self.points_grid);
            }
        }

        new_box
    }
}

impl<'a> PointsGrid {
    pub fn new(points: Vec<Vec<NotNan<f64>>>) -> Self {
        let d = points.len();
        let n = points[0].len();
        let zero = NotNan::new(0.0).unwrap();
        let one = NotNan::new(1.0).unwrap();
        let points: Vec<Vec<NotNan<f64>>> = points.iter().map(|points| {
            let zero = once(zero);
            let one = once(one);
            let points: Vec<_> = zero.chain(points.iter().copied()).chain(one).collect();
            points
        }).collect();
        let (ordered_points, reverse_order, order): (Vec<_>, Vec<Vec<PointIndex>>, Vec<Vec<PointOrder>>) = points.clone().into_iter().map(|points| {
            let indexed_points: Vec<_> = points.into_iter().enumerate().collect();
            let (points, reverse_order, order) = dedup_and_get_coord(indexed_points);
            let reverse_order = reverse_order.into_iter().map(|x| PointIndex(x)).collect();
            let order = order.into_iter().map(|x| PointOrder(x)).collect();
            (points, reverse_order, order)
        }).multiunzip();
        PointsGrid {
            n,
            d,
            points,
            ordered_points,
            order,
            reverse_order,
        }
    }

    pub fn one(&self) -> PointIndex {
        PointIndex(self.n + 1)
    }

    pub fn zero(&self) -> PointIndex {
        PointIndex(0)
    }

    pub fn point_order(&'a self, x: PointIndex) -> impl Iterator<Item = PointOrder> + 'a {
        self.order.iter().map(move |order| order[x.0])
    }

    pub fn new_origin_index(&'a self) -> Point<'a> {
        Point::new(
            repeat(self.zero()).take(self.d).collect(),
            self,
        )
    }

    pub fn new_max_index(&'a self) -> Point<'a> {
        Point::new(
            repeat(self.one()).take(self.d).collect(),
            self,
        )
    }

    pub fn get_permutation(&self, rng: &mut Rng) -> Vec<usize> {
        let mut order: Vec<_> = (0..self.d).collect();
        rng.shuffle(&mut order);
        order
    }

    pub fn generate_random_point(&'a self, rng: &mut Rng) -> RoughPoint<'a> {
        RoughPoint::new_random(&self, rng)
    }

    pub fn point_set(&self) -> impl Iterator<Item = PointIndex> {
        (1..self.n+1).map(PointIndex)
    }

    pub fn warp_with_pd_zo(&self, s: NotNan<f64>) -> NotNan<f64> {
        self.warp_with_pd(NotNan::new(0.0).unwrap(), NotNan::new(1.0).unwrap(), s)
    }

    pub fn warp_with_pd(&self, lower_bound: NotNan<f64>, upper_bound: NotNan<f64>, s: NotNan<f64>) -> NotNan<f64> {
        let d = self.d;
        let temp =
            NotNan::new(upper_bound.powi(d as i32) - lower_bound.powi(d as i32)).unwrap()
            *s
            + lower_bound.powi(d as i32);
        NotNan::new(temp.powf(1.0 / d as f64)).unwrap()
    }
}

// remove duplicate from indexed elements and return the permutation
pub fn dedup_and_get_coord<T: PartialEq>(stuff: Vec<(usize, T)>) -> (Vec<T>, Vec<usize>, Vec<usize>) {
    let mut reverse_order = Vec::<usize>::with_capacity(stuff.len());
    let mut order = Vec::<usize>::with_capacity(stuff.len());
    let mut just_things = Vec::<T>::with_capacity(stuff.len());
    stuff.into_iter().for_each(|(ii, vv)| {
        if let Some(v) = just_things.last() {
            if vv != *v {
                reverse_order.push(ii);
                just_things.push(vv);
            }
        } else {
            reverse_order.push(ii);
            just_things.push(vv);
        }
        order.push(reverse_order.len() - 1);
    });
    (just_things, reverse_order, order)
}

fn iterator_with_permutation<'b, T>(stuff: &'b [T], permutation: &'b [usize]) -> impl Iterator<Item = &'b T> + 'b {
    let iterator = permutation.iter().map(|i| &stuff[*i]);
    iterator
}

fn get_reveres_permutation(permutation: &[usize]) -> Vec<Option<usize>> {
    let mut reverse_permutation = vec![None; permutation.len()];
    for (i, j) in permutation.iter().enumerate() {
        reverse_permutation[*j] = Some(i);
    }
    reverse_permutation
}

fn reverse_iterator_with_permutation<T>(iterator: impl Iterator<Item = T>, reverse_permutation: &[usize]) -> Vec<Option<T>> {
    let mut ans = std::iter::repeat_with(|| Option::<T>::None).take(reverse_permutation.len()).collect::<Vec<_>>();
    for (i, v) in iterator.enumerate() {
        ans[reverse_permutation[i]] = Some(v);
    }
    ans
}

pub fn get_index_up_with_cmp<T, P>(search_vector: &[T], mut cmp: P) -> Option<usize>
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
    // ~~Because of the guard at the beginning, it's guaranteed that left + 1 is within bounds~~ I was wrong, if the length is one it will return out out bound (edge cases ughhhhhh).
    if left + 1 == search_vector.len() {
        return None;
    }
    Some(left + 1)
}

pub fn get_index_up<T: Ord>(search_vector: &[T], v: &T) -> Option<usize> {
    get_index_up_with_cmp(search_vector, |other| other.cmp(v))
}

pub fn get_index_down_with_cmp<T, P>(search_vector: &[T], mut cmp: P) -> Option<usize>
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

pub fn get_index_down<T: Ord>(search_vector: &[T], v: &T) -> Option<usize> {
    get_index_down_with_cmp(search_vector, |other| other.cmp(v))
}

fn subset_to_mask(subset: &[usize], n: usize) -> Vec<bool> {
    let mut mask = vec![false; n];
    for i in subset {
        mask[*i] = true;
    }
    mask
}
