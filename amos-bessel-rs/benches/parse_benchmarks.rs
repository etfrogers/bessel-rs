use std::env;
use std::fs;
use std::path::Path;

fn parse_point_estimate(content: &str) -> Option<f64> {
    let key = "\"point_estimate\":";
    if let Some(idx) = content.find(key) {
        let remainder = &content[idx + key.len()..];
        let end_idx = remainder.find([',', '}']).unwrap_or(remainder.len());
        let val_str = remainder[..end_idx].trim();
        val_str.parse::<f64>().ok()
    } else {
        None
    }
}

fn get_estimate(group: &str, func: &str, impl_name: &str, dataset: &str) -> Option<f64> {
    let path_str = format!(
        "target/criterion/{}/{} ({})/{}/new/estimates.json",
        group, func, impl_name, dataset
    );
    let path = Path::new(&path_str);
    if let Ok(content) = fs::read_to_string(path) {
        parse_point_estimate(&content)
    } else {
        None
    }
}

fn report_f32_vs_f64() {
    let functions = ["J", "Y", "I", "K", "H1", "H2", "Ai", "Bi"];
    let datasets = ["Scattered", "Grid"];
    let group = "f32 vs f64 Performance";

    println!("Performance Comparison: f32 vs f64");
    println!("{:-<68}", "");
    println!(
        "{:<10} | {:<12} | {:<10} | {:<10} | {:<7}",
        "Function", "Dataset", "f32 (ns)", "f64 (ns)", "Speedup"
    );
    println!("{:-<68}", "");

    let mut total_speedup = 1.0;
    let mut count = 0;

    for func in &functions {
        for dataset in &datasets {
            if let (Some(t_f32), Some(t_f64)) = (
                get_estimate(group, func, "f32", dataset),
                get_estimate(group, func, "f64", dataset),
            ) {
                let speedup = t_f64 / t_f32;
                total_speedup *= speedup;
                count += 1;
                println!(
                    "{:<10} | {:<12} | {:<10.2} | {:<10.2} | {:.2}x",
                    func, dataset, t_f32, t_f64, speedup
                );
            }
        }
    }

    if count > 0 {
        println!("{:-<68}", "");
        let geom_mean = total_speedup.powf(1.0 / count as f64);
        println!(
            "Overall Geometric Mean Speedup (f32 over f64): {:.2}x\n",
            geom_mean
        );
    } else {
        println!(
            "No benchmark data found. Please run `cargo bench --bench precision_bench` first.\n"
        );
    }
}

fn report_rust_vs_fortran() {
    let functions = ["J", "Y", "I", "K"];
    let datasets = ["Scattered", "Grid"];
    let group = "Rust vs Fortran";

    println!("Performance Comparison: Rust vs Fortran");
    println!("{:-<68}", "");
    println!(
        "{:<10} | {:<12} | {:<10} | {:<14} | {:<7}",
        "Function", "Dataset", "Rust (ns)", "Fortran (ns)", "Speedup"
    );
    println!("{:-<68}", "");

    let mut total_speedup = 1.0;
    let mut count = 0;

    for func in &functions {
        for dataset in &datasets {
            if let (Some(t_rust), Some(t_fortran)) = (
                get_estimate(group, func, "Rust", dataset),
                get_estimate(group, func, "Fortran", dataset),
            ) {
                let speedup = t_fortran / t_rust;
                total_speedup *= speedup;
                count += 1;
                println!(
                    "{:<10} | {:<12} | {:<10.2} | {:<14.2} | {:.2}x",
                    func, dataset, t_rust, t_fortran, speedup
                );
            }
        }
    }

    if count > 0 {
        println!("{:-<68}", "");
        let geom_mean = total_speedup.powf(1.0 / count as f64);
        println!(
            "Overall Geometric Mean Speedup (Rust over Fortran): {:.2}x\n",
            geom_mean
        );
    } else {
        println!("No benchmark data found. Please run `cargo bench --bench bessel_bench` first.\n");
    }
}

fn main() {
    let args: Vec<String> = env::args().collect();
    let mut run_f32_f64 = true;
    let mut run_rust_fortran = true;

    if args.len() > 1 {
        run_f32_f64 = false;
        run_rust_fortran = false;
        for arg in &args[1..] {
            match arg.as_str() {
                "--f32-vs-f64" => run_f32_f64 = true,
                "--rust-vs-fortran" => run_rust_fortran = true,
                _ => {
                    println!("Unknown argument: {}", arg);
                    println!("Usage: ./parse_benchmarks [--f32-vs-f64] [--rust-vs-fortran]");
                    return;
                }
            }
        }
    }

    if run_rust_fortran {
        report_rust_vs_fortran();
    }

    if run_f32_f64 {
        report_f32_vs_f64();
    }
}
