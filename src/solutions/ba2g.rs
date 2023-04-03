/*
    Implements GibbsSampler.
    https://rosalind.info/problems/ba2g/
 */
use rand::{
    Rng,
    distributions::WeightedIndex,
};

static BASES: &'static str = "ACGT";

/// Returns the best motifs initialized to random k-mers of each DNA string.
/// 
/// * `dna` - A vector of strings, each representing a DNA sequence.
/// * `k` - A positive integer representing the length of a single motif sequence.
fn seed_best_motifs(dna: &[String], k: usize) -> Vec<String> { // TODO
    let best_motifs: Vec<String> = dna.iter()
        .map(|seq| {
            let r = rand::thread_rng().gen_range(0..=seq.len() - k);
            seq[r..r+k].to_string()
        }).collect();

    best_motifs
}

/// Returns the profile matrix of a motif set according to Laplace's rule of succession.
/// 
/// * `motifs` - A vector of strings, each representing a motif.
fn gen_profile(motifs: &[String]) -> Vec<Vec<f32>> { // TODO
    let columns: Vec<Vec<char>> = (0..motifs[0].len())
        .map(|i| {
            motifs.iter()
                .map(|row| row.chars().nth(i).unwrap()).collect()
        }).collect();
    
    let profile: Vec<Vec<f32>> = columns.iter()
        .map(|col| {
            BASES.chars()
                .map(|base| (col.iter()
                    .filter(|&&c| c == base)
                    .count() as f32 + 1.0) / (col.len() as f32 + 4.0)).collect()
        }).collect();

    profile
}

/// Returns the profile-randomly generated k-mer in a sequence.
/// 
/// * `seq` - A single DNA string.
/// * `k` - A positive integer representing the length of a single motif sequence.
/// * `profile` - A vector of vectors of floats representing the profile matrix of a motif set.
fn profile_random_kmer(seq: &str, k: usize, profile: &[Vec<f32>]) -> String {
    let num_kmers = seq.len() - k + 1;
    let mut kmer_probs = vec![1.0; num_kmers];

    for start in 0..num_kmers {
        let kmer = &seq[start..start+k];
        for (col, base) in kmer.chars().enumerate() {
            kmer_probs[start] *= profile[col][BASES.find(base).unwrap()];
        }
    }

    let i = rand::thread_rng().sample(&WeightedIndex::new(&kmer_probs).unwrap());

    seq[i..i+k].to_string()
}

/// Returns the consensus string of a motif set.
/// 
/// * `motifs` - A vector of strings, each representing a motif.
fn gen_consensus(motifs: &[String]) -> String { // TODO
    let columns: Vec<Vec<char>> = (0..motifs[0].len())
        .map(|i| {motifs.iter()
            .map(|row| row.chars().nth(i).unwrap()).collect()
        }).collect();
        
    let col_counts: Vec<Vec<i32>> = columns.iter()
        .map(|col| {BASES.chars()
            .map(|base| col.iter().filter(|&&c| c == base).count() as i32).collect()
        }).collect();

    let consensus: String = col_counts.iter()
        .map(|counts| BASES.chars()
            .nth(counts.iter().enumerate()
                .max_by_key(|&(_, count)| count).unwrap().0).unwrap()).collect();
        
    consensus
}

/// Returns the hamming distance between a motif and a consensus string.
/// 
/// * `motif` - A string representing a motif.
/// * `consensus` - A string representing a consensus string.
fn hamming_distance(motif: &str, consensus: &str) -> i32 {
    let mut hdist = 0;

    for (m, c) in motif.chars().zip(consensus.chars()) {
        if m != c {
            hdist += 1;
        }
    }

    hdist
}

/// Returns the score of a motif set against a consensus string.
/// 
/// * `motifs` - A vector of strings, each representing a motif.
/// * `consensus` - A string representing a consensus string.
fn score(motifs: &[String], consensus: &str) -> i32 { // TODO
    let mut score = 0;

    for motif in motifs.iter() {
        score += hamming_distance(motif, consensus);
    }

    score
}

/// Returns the best motifs and their score by running GibbsSampler.
/// 
/// * `dna` - A vector of strings, each representing a DNA sequence.
/// * `k` - A positive integer representing the length of a single motif sequence.
/// * `t` - A positive integer representing the number of DNA sequences.
/// * `r` - A positive integer representing the number of iterations in the Gibbs sampling procedure.
pub fn gibbs_sampler(dna: Vec<String>, k: usize, t: usize, r: usize) -> (Vec<String>, i32) {
    let mut best_motifs = seed_best_motifs(&dna, k);
    let mut motifs = best_motifs.clone();
    let mut best_motifs_score = score(&best_motifs, &gen_consensus(&best_motifs));

    for _j in 0..r {
        let i = rand::thread_rng().gen_range(0..=t-1);
        motifs.remove(i);
        motifs.insert(i, profile_random_kmer(&dna[i], k, &gen_profile(&motifs)));

        let motifs_score = score(&motifs, &gen_consensus(&motifs));
        if motifs_score < best_motifs_score {
            best_motifs = motifs.clone();
            best_motifs_score = motifs_score;
        }     
    }

    (best_motifs, best_motifs_score)
}
