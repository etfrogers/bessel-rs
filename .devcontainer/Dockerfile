FROM ubuntu:22.04

ENV DEBIAN_FRONTEND=noninteractive

# Install build tools, compilers, and debuggers
RUN apt update && apt install -y \
    build-essential \
    cmake \
    gdb \
    clang \
    lldb \
    git \
    gfortran \
    curl \
    pipx \
    && rm -rf /var/lib/apt/lists/*

RUN curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs | bash -s -- -y --no-modify-path
RUN pipx install fortls
# RUN export PATH="/root/.cargo/bin:${PATH}"
# Add .cargo/bin to PATH
ENV PATH="/root/.cargo/bin:/root/.local/bin:${PATH}"
WORKDIR /workspace
