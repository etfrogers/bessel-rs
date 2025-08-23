cargo llvm-cov --html "$@"
cargo llvm-cov report --lcov --output-path lcov.info
open target/llvm-cov/html/index.html