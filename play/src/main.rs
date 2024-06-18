use itertools::izip;

fn main() {
    let mut yp_sn: Vec<_> = vec![0;5];
    let yp = vec![0;5];
    let x = vec![0;5];
    for (x, yp_sn, _yp) in izip!(x, yp_sn.iter_mut(), yp) {
        *yp_sn = x;
    }
}
