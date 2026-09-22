fn main() {
    println!("cargo::rustc-check-cfg=cfg(bessel_zeros_ref)");
    println!("cargo:rerun-if-changed=../bessel-zeros/Cargo.toml");
    if let Ok(content) = std::fs::read_to_string("../bessel-zeros/Cargo.toml")
        && content.contains("version = \"0.1.")
    {
        println!("cargo:rustc-cfg=bessel_zeros_ref");
    }
}
