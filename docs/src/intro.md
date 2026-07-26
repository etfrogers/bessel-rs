# Introduction

Welcome to the **Bessel-rs Performance & Accuracy Guide**.

`bessel-rs` is an idiomatic, pure-Rust translation of Amos' renowned complex Bessel function algorithms (originally written in Fortran). 

While porting legacy scientific code to memory-safe languages provides obvious security and ecosystem benefits, users of scientific libraries are rightfully concerned with two primary metrics before adopting a new crate:
1. **Performance**: Is it as fast as the legacy C/Fortran code? 
2. **Accuracy**: Does it preserve the mathematical precision and edge-case handling of the original algorithms?

This book serves to transparently document how `bessel-rs` performs across both of these domains.
