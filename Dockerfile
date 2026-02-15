FROM us-docker.pkg.dev/gemini-code-dev/gemini-cli/sandbox:0.1.1

USER root

# Install curl and other dependencies
# RUN apt update && apt install -y \
#     build-essential \
#     cmake \
#     gdb \
#     clang \
#     lldb \
#     git \
#     gfortran \
#     curl \
#     pipx \
#     && rm -rf /var/lib/apt/lists/*
RUN apt-get update && apt-get install -y \
    curl \
    build-essential \
    gcc \
    git \
    gdb \
    cmake \
    clang \
    lldb \
    gfortran \
    && rm -rf /var/lib/apt/lists/*

ENV RUSTUP_HOME=/home/node/.rustup \
    CARGO_HOME=/home/node/.cargo \
    PATH="/home/node/.cargo/bin:$PATH"
# Create the folders and fix permissions BEFORE installing Rust
RUN mkdir -p $RUSTUP_HOME $CARGO_HOME && chown -R node:node /home/node

USER node

# Install rustup
RUN curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs | sh -s -- -y
# add rustup to the path so we can use it in the next steps
# RUN . "/home/node/.cargo/env"
RUN /home/node/.cargo/bin/rustup toolchain install nightly
RUN /home/node/.cargo/bin/rustup component add clippy rustfmt --toolchain stable && \
    /home/node/.cargo/bin/rustup component add clippy rustfmt --toolchain nightly

# Add cargo to the path
# ENV PATH="/home/node/.cargo/bin:${PATH}"

# Set the working directory
WORKDIR /workspace


