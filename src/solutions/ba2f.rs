/*
    Implements RandomizedMotifSearch.
    https://rosalind.info/problems/ba2f/
 */
use rand::Rng;

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

/// Returns the profile-most probable k-mer in a sequence.
/// 
/// * `seq` - A single DNA string.
/// * `k` - A positive integer representing the length of a single motif sequence.
/// * `profile` - A vector of vectors of floats representing the profile matrix of a motif set.
fn profile_probable_kmer(seq: &str, k: usize, profile: &[Vec<f32>]) -> String {
    let mut max_prob = 0.0;
    let mut most_prob_kmer = &seq[0..k];

    for start in 0..seq.len()-k+1 {
        let kmer = &seq[start..start+k];
        let mut kmer_prob = 1.0;
        for (col, base) in kmer.chars().enumerate() {
            kmer_prob *= profile[col][BASES.find(base).unwrap()];
        }
        if kmer_prob > max_prob {
            max_prob = kmer_prob;
            most_prob_kmer = kmer;
        }
    }

    String::from(most_prob_kmer)
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

/// Returns the best motifs and their score by running RandomizedMotifSearch.
/// 
/// * `dna` - A vector of strings, each representing a DNA sequence.
/// * `k` - A positive integer representing the length of a single motif sequence.
/// * `t` - A positive integer representing the number of DNA sequences.
pub fn randomized_motif_search(dna: Vec<String>, k: usize, t: usize) -> (Vec<String>, i32) {
    let mut best_motifs = seed_best_motifs(&dna, k);
    let mut best_motifs_score = score(&best_motifs, &gen_consensus(&best_motifs));

    for _j in 0..1000 {
        let mut motifs = Vec::with_capacity(t);
        for i in 0..t {
            motifs.push(profile_probable_kmer(&dna[i], k, &gen_profile(&best_motifs)));
        }
        let motifs_score = score(&motifs, &gen_consensus(&motifs));
        if motifs_score < best_motifs_score {
            best_motifs = motifs;
            best_motifs_score = motifs_score;
        } else {
            break
        }
    }
    
    (best_motifs, best_motifs_score)
}