use ordered_float::NotNan;

// remove duplicate from indexed elements and return the permutation
pub fn dedup_and_get_coord<T: PartialEq>(mut stuff: Vec<(usize, T)>) -> (Vec<usize>, Vec<(usize, T)>) {
    let mut resulting_map: Vec<usize> = Vec::with_capacity(stuff.len());
    stuff.dedup_by(|a, b| {
        match (a.1 == b.1) {
            true => {
                resulting_map[b.0] = a.0;
                return true;
            },
            false => {
                resulting_map[b.0] = b.0;
                return false
            }
        }
    });
    return (resulting_map, stuff);

}

fn transpose<T>(v: Vec<Vec<T>>) -> Vec<Vec<T>> {
    assert!(!v.is_empty());
    let len = v[0].len();
    let mut iters: Vec<_> = v.into_iter().map(|n| n.into_iter()).collect();
    (0..len)
        .map(|_| {
            iters
                .iter_mut()
                .map(|n| n.next().unwrap())
                .collect::<Vec<T>>()
        })
        .collect()
}

pub fn process_coord_data(points: Vec<Vec<NotNan<f64>>>, d: usize, n: usize) -> (Vec<Vec<NotNan<f64>>, Vec<Vec<NotNan<f64>>) {
    let n_dimensions = d;
    let n_points = n;

    // let coordinate: Vec<usize> = (0..d).collect();
    
    let (coord, point_index): (Vec<Vec<NotNan<f64>>>, Vec<Vec<usize>>) = (0..n_dimensions).map(|i| -> (Vec<NotNan<f64>>, Vec<usize>) {

        let zero = NotNan::new(0.0f64).unwrap();
        let zero = (0usize, zero);
        let one = NotNan::new(1.0f64).unwrap();
        let one = (n_points+1, one);
        let mut temp_coord: Vec<(usize, NotNan<f64>)> = std::iter::once(zero)
            .chain((0..n_points).map(|j| (j+1, points[j][i])))
            .chain(std::iter::once(one))
            .collect();

        // TODO: figure out if sorting is needed
        temp_coord.sort_by(|a, b| a.1.cmp(&b.1));

        let (index_map, temp_coord) = dedup_and_get_coord(temp_coord);
        
        let ans = temp_coord.into_iter().map(|a| a.1).collect();
        
        return (ans, index_map);

    }).unzip();

    // FIXME: remove this transpose
    return (coord, transpose(point_index));

}