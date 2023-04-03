use std::{
    error::Error,
    path::Path,
    io::{BufReader, BufRead},
    fs::File,
};

mod solutions;

#[derive(Debug)]
pub struct Config {
    pub problem: String,
}

impl Config {
    pub fn build(
        mut args: impl Iterator<Item=String>
    ) -> Result<Config, &'static str> {
        args.next();
        let problem = match args.next() {
            Some(arg) => arg,
            None => return Err("Didn't get a file path"),
        };
        
        Ok(Config { problem })
    }
}

pub fn solve(problem: &str) {
    match problem {
        "RandomizedMotifSearch" => solve_ba2f(),
        "GibbsSampler" => solve_ba2g(),
        _ => panic!("Unknown function name: {}", problem),
    }
}

fn open(filename: &str) -> Vec<Result<String, std::io::Error>> {
    let path = Path::new(filename);
    let file = match File::open(&path) {
        Err(why) => panic!("Couldn't open {}: {}", path.display(), <dyn Error>::to_string(&why)), 
        Ok(file) => file,
    };
    let reader = BufReader::new(file);
    
    reader.lines().collect()
}

fn solve_ba2f() {
    println!("Solving RandomizedMotifSearch (ba2f)");
    let lines = open("inputs/ba2f.txt");
    let kt: Vec<&str> = if let Ok(s) = &lines[0] {
        s.split_whitespace().collect()
    } else {Vec::new()};
    let k = kt[0].parse::<usize>().unwrap();
    let t = kt[1].parse::<usize>().unwrap();
    let dna: Vec<String> = lines[1..lines.len()].into_iter()
        .map(|r| r.as_ref().unwrap().to_string()).collect();
    let (best_motifs, best_motifs_score) = solutions::ba2f::randomized_motif_search(dna, k, t);
    for m in best_motifs {
        println!("{m}");
    }
    println!("{best_motifs_score}");
}

fn solve_ba2g() {
    println!("Solving GibbsSampler (ba2g)");
    let lines = open("inputs/ba2g.txt");
    
}