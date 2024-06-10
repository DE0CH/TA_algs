pub mod common;
use common::hello;
fn main() {
    println!("Hello, world!");
    hello(String::from("Alice"));
}
