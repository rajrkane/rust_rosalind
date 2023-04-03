use std::{
    env, 
    process,
};
use rust_rosalind::Config;

fn main() {
    let config = Config::build(env::args()).unwrap_or_else(|err| {
        eprintln!("Problem parsing arguments: {}", err);
        process::exit(1);
    });

    rust_rosalind::solve(&config.problem);
}
