////use std::env;
use std::fs::File;
use std::io::{self, BufRead};
use std::path::Path;
///////////////////////
use needletail::parse_fastx_file;
use clap::Parser;
use std::convert::TryInto;
use substring::Substring;
use std::collections::HashMap;
use regex::Regex;
use distance::*;
////////////////////////////////

// Command line arguments -- configured in the structure name Cli
#[derive(Parser, Debug)]
#[command(version, about, long_about = None)]
struct Cli {
    /// This is a string containing barcode information
    #[arg(short, long, default_value="0_7+11_16+20_25+31_38")]
    scrna_barcode: Option<String>,
    /// The path to the FASTQ file to read to extract barcodes
    #[arg(short,long)]
    fastq: std::path::PathBuf,
    /// barcode whitelist file
    #[arg(short,long)]
    barcode_file:  std::path::PathBuf,
}

// from a vector of scores , return back overall score for candidate
fn score_candidate(list: &Vec<i32>) -> f32 {
    // create boolean vector
    let bool_iter = list.iter().map(|&i| i as f32 == 0.0).into_iter();
    // sum based on true / false
    let num_hits: f32 = bool_iter.map(|b| if b { 1.0 } else { 0.0 }).sum();
    let total_diff: f32 = list.iter().map(|&i| i as f32 ).sum();
    let mut score : f32 = 0.0;
    if list.len() as f32 == num_hits {
        score = num_hits;
    } else if num_hits == 0.0 {
        score = 0.00001/total_diff;
    } else {
        score = num_hits/total_diff;
    }
    return score;
}


// can variable phase 0-3 bases from beginning of read when trying to extract barcode
// def phased_check(original_sequence,barcode_string,lookup_list):
//    linker_sequences = ['ATG','CAG','TCGAG']
//    ### initialize scores
//    my_best_hit = dict()
//    dummy_barcode = "N" * len(original_sequence)
//    my_best_hit['barcode_seq'] = dummy_barcode
//    my_best_hit['best_barcode'] = dummy_barcode
//    my_best_hit['best_dists'] = [100000,100000,100000,100000]
//    my_best_hit['best_tiered_seqs'] = ["N" * 8, "N" * 6,"N" * 6,"N" * 8]
//    /// CURRRENTLY will not check for linker sequences, but will return barcode with the lowest hamming distance
//    phase_lag = [0 , 1, 2, 3]
//    keys_to_check = ['block1_dist','block2_dist','block3_dist','block4_dist']
//    for i in phase_lag:
//        candidate_barcodes = get_barcode_sequences(original_sequence, barcode_string,pad = i)
//        candidate_linkers = get_linker_sequences(original_sequence, barcode_string,pad = i)
//        linkers_match = candidate_linkers['linker_blocks'] == linker_sequences
//        candidate_dists = compute_hamming_by_tier(candidate_barcodes['barcode_blocks'],lookup_list)
//        candidate_scores = []
//        for k in keys_to_check:
//            candidate_scores.append(candidate_dists[k])
//        logging_statement(f"BEST_HIT:\t{my_best_hit['best_dists']}\tCANDIDATE_HIT:\t{candidate_scores}\t{i}")  
//        if score_candidate(candidate_scores) > score_candidate(my_best_hit['best_dists']) or linkers_match is True:
//            #for idx,j in enumerate(keys_to_check):
//            my_best_hit['best_dists'] = candidate_scores
//            my_best_hit['best_tiered_seqs'] = [candidate_dists['block1_str'],candidate_dists['block2_str'],candidate_dists['block3_str'],candidate_dists['block4_str']]
//            my_best_hit['barcode_seq'] = candidate_barcodes['barcode_str']
//            my_best_hit['best_barcode'] = "".join( my_best_hit['best_tiered_seqs'])
//    return my_best_hit

 fn phased_check(read_sequence: String, barcode_string: String,whitelist_map: HashMap<String,Vec<String>>,phase_values: Vec<i32>) -> HashMap<std::string::String, Vec<String>>{
    // create hash map to track best score
    let mut my_best_hit: HashMap<String, Vec<String>> = HashMap::new();
    let mut barcode_seq: String = "NNNNN".to_string();
    // my_best_hit.entry("barcode_seq").or_insert(BarcodeValue::String("NNNNN".to_string()));
    let mut best_barcode: String = "NNNNN".to_string();
    // my_best_hit.entry("best_barcode").or_insert(BarcodeValue::String("NNNNN".to_string()));
    let mut best_dists: Vec<i32> = vec![10000,10000,10000,10000];
    // my_best_hit.entry("best_dists").or_insert(BarcodeValue::VecofNums([10000,10000,10000,10000]));
    let mut best_tiered_seqs: Vec<String> = vec!["NNNNN".to_string(),"NNNNN".to_string(),"NNNNN".to_string(),"NNNNN".to_string()];
    //my_best_hit.entry("best_tiered_seqs").or_insert(BarcodeValue::VecofString(["NNNNN","NNNNN","NNNNN","NNNNN"]));
    let mut keys_to_check: Vec<_>  = whitelist_map.keys().collect();
    keys_to_check.sort();
    // outer loop on all phase values
    for index in 0..phase_values.len(){
        let my_barcode_candidate = get_barcode_str(&read_sequence,&barcode_string,phase_values[index]);
        let mut barcode_strings: Vec<String> = Vec::new();
        let mut barcode_dists: Vec<i32> = Vec::new();
        let mut all_barcode_strings: Vec<String> = Vec::new();
        for k in 0..keys_to_check.len(){
        // second loop on all barcode blocks
            let mut my_block = whitelist_map.get(keys_to_check[k]).unwrap();
            let mut barcode_to_check: String = my_barcode_candidate.get(keys_to_check[k]).unwrap().to_string();
            all_barcode_strings.push(barcode_to_check.clone());
            // println!("barcode_to_check: {:?}",barcode_to_check);
            // println!("my_block[0]: {:?}",my_block[0]);
            let my_block_match = get_best_match(barcode_to_check,my_block);
            let my_block_match_keys: Vec<_>  = my_block_match.keys().collect();
            let my_dist = my_block_match.get(my_block_match_keys[0]).unwrap();
            barcode_strings.push(my_block_match_keys[0].to_string());
            barcode_dists.push(*my_dist);

            // println!("my_block_match_keys: {:?}",my_block_match_keys[0]);
            // println!("my_block_match_keys_hamming {:?}",my_block_match.get(my_block_match_keys[0]).unwrap());
        }
        // once we collect dists for all scores
        let current_score = score_candidate(&barcode_dists);
        let best_score = score_candidate(&best_dists);
        if current_score > best_score {
            best_dists = barcode_dists;
            //my_best_hit.insert("best_dists", barcode_dists);
            best_tiered_seqs = barcode_strings.clone();
            //my_best_hit.insert("best_tiered_seqs", barcode_strings);
            barcode_seq = all_barcode_strings.join("").to_string();
            //my_best_hit.insert("barcode_seq",all_barcode_strings.join("").to_string());
            best_barcode = barcode_strings.join("").to_string();
            //my_best_hit.insert("best_barcode",barcode_strings.join("")).to_string();
        }
    }
    let final_dists = best_dists.iter().map(|&i| i.to_string()).collect();
    my_best_hit.entry("barcode_seq".to_string()).or_insert_with(Vec::new).push(barcode_seq);
    my_best_hit.entry("best_tiered_seqs".to_string()).or_insert(best_tiered_seqs);
    my_best_hit.entry("best_barcode".to_string()).or_insert_with(Vec::new).push(best_barcode);
    my_best_hit.entry("best_dists".to_string()).or_insert(final_dists);
    return my_best_hit;
 }


/// Return the best hit(barcode) and distance score (hamming)
fn get_best_match(str1: String, str_array: &Vec<String>) ->  HashMap<std::string::String, i32> {
    let mut best_barcode: String = "NNNNNNNN".to_string();
    let mut best_distance = 100000;
    let mut bestmatch_map = HashMap::new();
    for index in 0..str_array.len(){
        let str2 = &str_array[index];
        // println!("str2: {str2}");
        let distance = hamming(&str1, &str2).unwrap();
        if  distance < best_distance { 
            best_distance = distance.to_owned();
            best_barcode =  str2.to_owned();
        }
    }

    bestmatch_map.entry(best_barcode.to_string())
    .or_insert(best_distance.to_owned() as i32);    
    return bestmatch_map;
}
// The output is wrapped in a Result to allow matching on errors.
// Returns an Iterator to the Reader of the lines of the file.
fn read_lines<P>(filename: P) -> io::Result<io::Lines<io::BufReader<File>>>
where P: AsRef<Path>, {
    let file = File::open(filename)?;
    Ok(io::BufReader::new(file).lines())
}

fn get_barcode_whitelist(whitelist_file:  &std::path::PathBuf) -> HashMap<std::string::String, Vec<std::string::String>> {
    let mut whitelist_map = HashMap::new();
    let _re = Regex::new(r"^#").unwrap();
    let mut key: String = "UNKNOWN".to_string();

    // File hosts must exist in current path before this produces output
    if let Ok(lines) = read_lines(whitelist_file) {
        // Consumes the iterator, returns an (Optional) String
        for line in lines {
            if let Ok(data) = line {
                let values: Vec<&str> = data.split('-').collect();
                match values.len() {
                  2 => {
                    key = values[values.len()-1].to_string();
                  },
                  1 => {
                    whitelist_map.entry(key.to_string())
                    .or_insert_with(Vec::new)
                    .push(data.to_owned());                    
                  },                 
                   _ => panic!("Invalid input line {}", data),
                };
            }
        }
    }
    return whitelist_map;
}

// function to get barcode string subset
fn get_barcode_str(og_str: &str,barcode_str: &String, padding: i32) -> HashMap<std::string::String, std::string::String> {
    let mut barcode_map = HashMap::new();
    // split barcode string and turn it into barcode string
    let parts = barcode_str.as_str().split("+");
    //let mut og_str_indices = og_str.char_indices();
    let collection = parts.clone().collect::<Vec<&str>>();
    let mut index_num: i32 = 0;
    for part in collection {
        index_num += 1;
        let my_key = "Block".to_owned() + &index_num.to_string();
        let interval_split = part.split("_");
        let collection1 = interval_split.clone().collect::<Vec<&str>>();
        let numbers: Result<Vec<i32>, _> = collection1.iter().map(|x| x.parse()).collect();
        //dbg!(&numbers);
        let indexes = numbers.unwrap();
        let start_index: i32 = indexes[0] + padding;
        // println!("start_index: {start_index}");
        let end_index: i32 = indexes[1] + padding + 1;
        // println!("end_index: {end_index}");
        let substring = og_str.substring(start_index.try_into().unwrap(),end_index.try_into().unwrap());
        // println!("{:?}",substring);
        barcode_map.entry(my_key.to_string())
        .or_insert(substring.to_owned());  
    }
    return barcode_map;
}

fn main() {
    let args = Cli::parse();
    /////////////////////////////////////
    let my_str = "CCCAAAAAAAAAAAAAAAAAAABBBBBBBBBBBBBBBBBBBZZZ";
    let my_barcode_str_clone = args.scrna_barcode.clone().unwrap();
    let my_barcode_candidate = get_barcode_str(&my_str,&args.scrna_barcode.unwrap(),1);
    let my_whitelist = get_barcode_whitelist(&args.barcode_file);
    let mut iter_keys: Vec<_>  = my_whitelist.keys().collect();
    iter_keys.sort();
    // println!("iter_keys: {:?}",iter_keys);

    //let first_block= my_whitelist.get(iter_keys[0]).unwrap();
    //let barcode_to_check: String = my_barcode_candidate.get(iter_keys[0]).unwrap().to_string();
    // println!("barcode_to_check: {:?}",barcode_to_check);

    //let my_first_block_match = get_best_match(barcode_to_check,first_block);
    //let my_first_block_match_keys: Vec<_>  = my_first_block_match.keys().collect();
    // println!("my_first_block_match_keys: {:?}",my_first_block_match_keys[0]);
    // println!("my_first_block_match_keys_hamming {:?}",my_first_block_match.get(my_first_block_match_keys[0]).unwrap());

    //let test_scores: Vec<i32> = vec![1,1,1,1];
    //let my_score = score_candidate(&test_scores);
    // println!("test_scores: {:?} my_score {:?}",test_scores,my_score);

    let my_phase_values: Vec<i32> = vec![0,1,2,3];
    let my_final_check = phased_check(my_str.to_string(), my_barcode_str_clone,my_whitelist,my_phase_values);
    //  println!("my_final_check: {:?}",my_final_check);
    
    let best_tiered_seqs = my_final_check.get("best_tiered_seqs").unwrap();
    let barcode_seq = my_final_check.get("barcode_seq").unwrap();
    let best_barcode = my_final_check.get("best_barcode").unwrap();
    let best_dists = my_final_check.get("best_dists").unwrap();
    // Use map() to parse each string to an i32
    let int_vec: Result<Vec<i32>, _> = best_dists.iter().map(|s| s.parse::<i32>()).collect();
    let total_dist: i32 = int_vec.unwrap().iter().sum();
    println!("{},{},{},{},{},{},{},{},{},{},{}",barcode_seq[0],best_barcode[0],best_dists[0],best_tiered_seqs[0],best_dists[1],best_tiered_seqs[1],best_dists[2],best_tiered_seqs[2],best_dists[3],best_tiered_seqs[3],total_dist);
    std::process::exit(0);
    //////////////////////
    let mut n: u32 = 0;
    let mut slen: u64 = 0;
    let mut qlen: u64 = 0;
    let mut reader = parse_fastx_file(&args.fastq).expect("valid path/file");
    let header_line_str = "read_name,barcode_sequence_extracted,closest_whitelist_barcode,block1_dist,block1_str,block2_dist,block2_str,block3_dist,block3_str,block4_dist,block4_dist,total_dist";
    println!("{}",header_line_str);
    while let Some(record) = reader.next() {
        let seqrec = record.expect("invalid record");
        n += 1;
        slen += seqrec.seq().len() as u64;
        if let Some(qual) = seqrec.qual() {
            qlen += qual.len() as u64;
        }
    }
    //println!("{}\t{}\t{} found in fastq {}", n, slen, qlen,args.fastq.display());
    //println!("scrna_barcode: {:?}, fastq: {:?}", args.scrna_barcode, args.fastq);
    std::process::exit(0)
}