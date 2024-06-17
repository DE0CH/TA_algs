use clap::{arg, value_parser, Command};
use ordered_float::NotNan;
use std::fs::File;
use std::io::{BufRead, BufReader};

fn read_input(input_file: &str, d: usize, n: usize) -> Vec<Vec<NotNan<f64>>> {
    let file = File::open(input_file).unwrap();
    let mut file = BufReader::new(file);
    let ans = (0..n).map(|i| {
        let mut input = String::new();
        file.read_line(&mut input).expect(&format!("Failed to read line {}. Perhaps there aren't enough lines for {} points", i + 1, n));
        let point: Vec<_> = input.split_whitespace().map(|x| {
            x.parse::<f64>()
            .expect(&format!("Failed parse floating point on line {}", i + 1))
        }).map(|x| {
            NotNan::new(x)
            .expect(&format!("Floating point is NaN on line {}", i + 1))
        }).collect();
        if point.len() != d {
            panic!("Point line {} does not have the correct dimension", i + 1);
        }
        point
    }).collect();
    let mut input: String = String::new();
    if let Ok(_) = file.read_line(&mut input) {
        for x in input.split_whitespace() {
            x.parse::<f64>()
            .expect_err(&format!("There appears to be more numbers in the file after n (={}) lines", n));
        }
    }
    ans
}

fn check_points(raw_points: &Vec<Vec<NotNan<f64>>>) {
    raw_points.iter().enumerate().for_each(|(i, p)| {
        p.iter().for_each(|x| {
            let one = NotNan::new(1.0).unwrap();
            let zero = NotNan::new(0.0).unwrap();
            if *x > one || *x < zero {
                panic!("The {}'th point is not in the unit cube", i + 1);
            }
        })
    });
}

pub fn get_raw_points() -> (Vec<Vec<NotNan<f64>>>, u64, u64) {
    let matches = Command::new("")
        .arg(
            arg!(iterations: -i --iterations [iterations] "Number of Iterations")
            .value_parser(value_parser!(u64))
            .default_value("100000")
        ).arg(
            arg!(seed: -s --seed [seed] "Seed for the rng")
            .value_parser(value_parser!(u64))
            .default_value("42")
        ).arg(
            arg!(dimension: <dimension> "Dimension")
            .value_parser(value_parser!(usize))
        ).arg(
            arg!(n_points: <n_points> "Number of points")
            .value_parser(value_parser!(usize))
        ).arg(
            arg!(input: <input> "Input file")
        ).get_matches();
    let iterations = *matches.get_one::<u64>("iterations").unwrap();
    let d = *matches.get_one::<usize>("dimension").unwrap();
    let n = *matches.get_one::<usize>("n_points").unwrap();
    let seed = *matches.get_one::<u64>("seed").unwrap();
    let input_file = matches.get_one::<String>("input").unwrap();
    println!("Number of iterations: {}", iterations);
    println!("Dimension: {}", d);
    println!("Number of points: {}", n);
    let raw_points = read_input(input_file, d, n);
    check_points(&raw_points);
    let mut processed_points = Vec::<Vec<NotNan<f64>>>::with_capacity(n);
    (0..d).for_each(|id| {
        let mut points = Vec::<NotNan<f64>>::with_capacity(n);
        raw_points.iter().for_each(|p| {
            points.push(p[id]);
        });
        processed_points.push(points);
    });
    (processed_points, seed, iterations)
}
