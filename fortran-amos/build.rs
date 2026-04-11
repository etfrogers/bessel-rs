fn main() {
    // 1. Tell Cargo to rerun this script if the Fortran files change
    println!("cargo:rerun-if-changed=src/amos_iso_c_fortran_wrapper.f90");
    println!("cargo:rerun-if-changed=src/machine.for");
    println!("cargo:rerun-if-changed=src/zbesh.for");

    // 2. Configure the Fortran compiler
    cc::Build::new()
        .compiler("gfortran")
        .files([
            "src/amos_iso_c_fortran_wrapper.f90",
            "src/machine.for",
            "src/zbesh.for",
        ])
        // Basic flags
        .flag("-O0")
        .flag("-Wno-maybe-uninitialized")
        .flag("-Wno-compare-reals")
        .flag("-Wno-intrinsic-shadow")
        .flag("-Wno-do-subscript")
        .flag("-Wno-unused-dummy-argument")
        // Optional: Adding these might help stabilize the "NaN in CI" issue
        // by preventing aggressive floating-point reordering.
        .flag("-fno-unsafe-math-optimizations")
        .flag("-frounding-math")
        .flag("-fsignaling-nans")
        .compile("amos_testing"); // This creates libamos_testing.a
}
