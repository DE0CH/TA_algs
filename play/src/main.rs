#![feature(layout_for_ptr)]
use std::{alloc::Layout, ops::Drop};

trait Help {}

#[derive(Debug)]
struct Me {
    data: Vec<String>
}

impl Drop for Me {
    fn drop(&mut self) {
        println!("Dropping Me");
    }
}

impl Help for Me {}

fn main() {
    let me = Me { data: vec!["please".to_string(), "help".to_string()] };
    let me: Box<dyn Help> = Box::new(me);
    let p = Box::into_raw(me) as *mut Me;
    let inner = unsafe {
        std::mem::transmute::<Box::<dyn Help>, Box::<Me>>(me)
    };
    dbg!(&inner.data);
    assert_eq!(inner.data, vec!["hello", "world"]);
}
