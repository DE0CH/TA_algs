fn main() {
    let mut v: Vec<i32> = Vec::with_capacity(10);
    unsafe {v.set_len(10)};
    v[9] = 10;
    println!("{:?}", v);
}