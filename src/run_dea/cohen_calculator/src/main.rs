use csv::{Reader, ReaderBuilder, StringRecord};
use indicatif::ProgressBar;
use serde::{Deserialize, Serialize};
use std::collections::HashMap;
use std::fs::File;
use std::io::BufReader;
use std::path::PathBuf;
use std::{error::Error, iter::Iterator};
use clap::{Parser, command};

#[derive(Parser,Debug)]
#[command(author = "Luca Visentin", about= "Calculate cohen's D of expression values.")]
struct Args {
    /// Path to the expression matrix with the expression values
    expression_matrix: PathBuf,
    /// Path to the metadata file with the sample descriptors
    metadata_file: PathBuf,
    /// Path to the json file with the case/control metadata matches
    matches_json: PathBuf,
    /// Path and filename of the output file
    output_path: PathBuf
}

// TODO: Combine these into one with generics? But I don't know exactly how...
#[derive(Debug)]
struct Mapper {
    case: Vec<String>,
    control: Vec<String>,
}

#[derive(Debug)]
struct Indexer {
    case: Vec<usize>,
    control: Vec<usize>,
}

#[derive(Serialize, Deserialize, Clone, Debug)]
struct TcgaToGtex {
    case_variable: String,
    control_variable: String,
    case_value: String,
    control_value: String,
}

fn main() {
    // Parse the command-line arguments.
    let args: Args = Args::parse();

    let samples: StringRecord = csv::ReaderBuilder::new().delimiter(b'\t').from_path(&args.expression_matrix).expect("Couldn't read header").headers().expect("expected a header").to_owned();

    println!("Reading Matches file...");
    let matches: HashMap<String, TcgaToGtex> = serde_json::from_reader(BufReader::new(
        File::open(&args.matches_json).expect("Cannot find matches file"),
    ))
    .expect("Invalid Json?");

    // Make a container for the new matches.
    let mut boxes: HashMap<String, Mapper> = HashMap::new();
    // Cloning so i can use into_iter() - unsure on how it works with .iter()
    // alone... So I don't do it.
    let match_iter = matches.clone().into_iter();
    for (match_label, _) in match_iter {
        boxes.insert(
            match_label,
            Mapper {
                case: Vec::new(),
                control: Vec::new(),
            },
        );
    }

    // Find the case and control samples
    println!("Generating matches...");
    let mut meta_reader = Reader::from_path(&args.metadata_file).unwrap();
    let pb = ProgressBar::new(18000); // some arbitrary len
    // This deserializes a record based on the header to a HashMap of String: String
    // With the header: Value key. This is to make it easier to grow the matches
    // 'boxes'.
    for line in meta_reader.deserialize() {
        let record: HashMap<String, String> = line.expect("Failed to read line");

        let id = record["Sample"].to_owned();
        if (id == "") | (! samples.clone().iter().any(|x| x == id)){
            // Skip if this has no IDs
            continue;
        }

        let match_iter = matches.clone().into_iter();
        for (match_label, item) in match_iter {
            // If the ID has the variable/value that we want, push it to the
            // vector of the correct box.

            // Case ----
            if record[&item.case_variable] == item.case_value {
                boxes
                    .get_mut(&match_label)
                    .unwrap()
                    .case
                    .push(id.clone());
            }

            // Control ----
            if record[&item.control_variable] == item.control_value {
                boxes
                    .get_mut(&match_label)
                    .unwrap()
                    .control
                    .push(id.clone());
            }
        }
        pb.inc(1);
    }

    println!("Computing cohen's D...");
    match process_csv(&args.expression_matrix, boxes, &args.output_path) {
        Ok(_) => println!("Done!"),
        Err(e) => println!("Something went wrong: {:?}", e),
    }
}

/// TODO: The mean, var and cohen functions would be more useful if we made them
/// generic. But I don't know how to do that (very well), so I leave them be.

/// Calculate the mean of the values in the input vector
fn mean(data: Vec<f32>) -> Option<f32> {
    let sum = data.iter().sum::<f32>();
    let count = data.len();

    match count {
        positive if positive > 0 => Some(sum / count as f32),
        _ => None,
    }
}

/// Calculate the variance of the values in the input vector
fn var(data: Vec<f32>) -> f32 {
    let data_mean = mean(data.clone()).unwrap();
    let count = &data.len();

    let variance = data
        .iter()
        .map(|value| {
            let diff = data_mean - (*value as f32);

            diff * diff
        })
        .sum::<f32>()
        / *count as f32;

    variance
}

/// Calculate cohen's D statistic from a case and control numeric vectors.
fn cohen(case: Vec<f32>, control: Vec<f32>) -> f32 {
    let n_case = case.len() as i32;
    let n_control = case.len() as i32;

    let pooled_var = (((n_case - 1) as f32 * var(case.clone())
        + (n_control - 1) as f32 * var(control.clone()))
        / (n_case + n_control - 2) as f32)
        .powf(0.5);

    if pooled_var == 0.0 {
        return 0.0;
    }

    (mean(case).unwrap() - mean(control).unwrap()) / pooled_var
}


fn process_csv(path: &PathBuf, boxes: HashMap<String, Mapper>, output_path: &PathBuf) -> Result<(), Box<dyn Error>> {
    let mut reader = ReaderBuilder::new().delimiter(b'\t').from_path(path)?;

    // We need the header for later, and we can also print the number of cols
    let header = reader.headers()?;
    println!("Read header. Found {} items", header.len());

    // We need to move from a Mapper to a Indexer. This does that, by finding
    // the indexes of the IDs in the mappers.
    let mut indexes: HashMap<String, Indexer> = HashMap::new();
    let biter = boxes.into_iter();

    // We iterate over the header values, saving the index of the found value.
    for (label, map) in biter {
        let case_labels = map
            .case
            .iter()
            .map(|x| {
                header
                    .iter()
                    .position(|i| i == *x)
                    .expect(format!("Missing case label {}", x).as_str())
            })
            .collect();

        let control_labels: Vec<usize> = map
            .control
            .iter()
            .map(|x| {
                header
                    .iter()
                    .position(|i| i == *x)
                    .expect(format!("Missing control label {}", x).as_str())
            })
            .collect();
        
        // Save the found indexes in the new "Indexer" structs
        indexes.insert(
            label,
            Indexer {
                case: case_labels,
                control: control_labels,
            },
        );
    }

    // We can now start calculating the statistic.

    // Initialize the output file
    let mut writer = csv::Writer::from_path(output_path)?;
    // Write the csv header first, just once. The first item must be "gene_id"
    // since we will write those first in every row too.
    let mut out_header = vec!["gene_id".to_string()];

    // We need to keep the order of the indexes static. So I clone them here by
    // collecting to a vector.
    let header_indexes: Vec<String> = indexes.keys().map(|x| x.to_string()).collect();
    out_header.extend(header_indexes.clone());
    writer.write_record(out_header)?;

    // All the PBs have arbitrary values, just to give them some space.
    let pb = ProgressBar::new(60000);
    for line in reader.records() {
        let record = line?;

        // Create a container for the result.
        let mut computed_cohen: HashMap<String, f32> = HashMap::new();
        // For every line, we need to iterate and compute the statistic.
        // here, "label" is the gene_id.
        for (label, dex) in indexes.iter() {
            // We know how many values we will find. So we can initialize the
            // number vectors with the correct size.
            let mut case: Vec<f32> = Vec::with_capacity(dex.case.len());
            let mut control: Vec<f32> = Vec::with_capacity(dex.control.len());

            // We have to parse the str from the csv reader to an f32 of the
            // values that we need.
            for x in &dex.case {
                case.push(record[*x].parse::<f32>().unwrap())
            }
            for x in &dex.control {
                control.push(record[*x].parse::<f32>().unwrap())
            }

            // Compute the statistic, and add it to the output hashmap.
            computed_cohen.insert(label.to_owned(), cohen(case, control));
        }

        // Create a container to build the output line that we need to write
        let mut result: Vec<String> = Vec::new();

        result.push(record[0].to_owned()); // push the gene id first
        // Push the values in the same order as the header that we saved first.
        result.extend(header_indexes.iter().map(|x| computed_cohen[x].to_string()));

        writer.write_record(result)?;

        pb.inc(1);
    }

    writer.flush()?;

    Ok(())
}
