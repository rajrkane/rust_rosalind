use std::{
    fs::{self, File}, 
    io::{BufReader, BufRead},
    path::Path,
    error::Error,
};

mod solutions;

fn main() {
    let path = Path::new("inputs/ba2f.txt");
    let file = match File::open(&path) {
        Err(why) => panic!("couldn't open {}: {}", path.display(), <dyn Error>::to_string(&why)), 
        Ok(file) => file,
    };
    let reader = BufReader::new(file);
    let lines: Vec<_> = reader.lines().collect();
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
