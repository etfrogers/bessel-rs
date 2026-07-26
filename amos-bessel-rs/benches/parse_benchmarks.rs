use std::fs;
use std::path::Path;

fn parse_point_estimate(content: &str) -> Option<f64> {
    let key = "\"point_estimate\":";
    if let Some(idx) = content.find(key) {
        let remainder = &content[idx + key.len()..];
        let end_idx = remainder.find(|c: char| c == ',' || c == '}').unwrap_or(remainder.len());
        let val_str = remainder[..end_idx].trim();
        val_str.parse::<f64>().ok()
    } else {
        None
    }
}

fn get_estimate(func: &str, prec: &str, dataset: &str) -> Option<f64> {
    let path_str = format!(
        "../target/criterion/f32 vs f64 Performance/{} ({})/{}/new/estimates.json",
        func, prec, dataset
    );
    let path = Path::new(&path_str);
    if let Ok(content) = fs::read_to_string(path) {
        parse_point_estimate(&content)
    } else {
        None
    }
}

fn main() {
    let functions = ["J", "Y", "I", "K", "H1", "H2", "Ai", "Bi"];
    let datasets = ["Scattered", "Grid"];

    println!("Performance Comparison: f32 vs f64");
    println!("{:-<60}", "");
    println!("{:<10} | {:<12} | {:<10} | {:<10} | {}", "Function", "Dataset", "f32 (ns)", "f64 (ns)", "Speedup");
    println!("{:-<60}", "");

    let mut total_speedup = 1.0;
    let mut count = 0;

    for func in &functions {
        for dataset in &datasets {
            if let (Some(t_f32), Some(t_f64)) = (
                get_estimate(func, "f32", dataset),
                get_estimate(func, "f64", dataset)
            ) {
                let speedup = t_f64 / t_f32;
                total_speedup *= speedup;
                count += 1;
                println!("{:<10} | {:<12} | {:<10.2} | {:<10.2} | {:.2}x", func, dataset, t_f32, t_f64, speedup);
            }
        }
    }

    if count > 0 {
        println!("{:-<60}", "");
        let geom_mean = total_speedup.powf(1.0 / count as f64);
        println!("Overall Geometric Mean Speedup (f32 over f64): {:.2}x", geom_mean);
    } else {
        println!("No benchmark data found. Please run `cargo bench --bench precision_bench` first.");
    }
}
