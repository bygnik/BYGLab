FROM ubuntu:22.04

RUN apt-get update && apt-get install -y --no-install-recommends \
    build-essential \
    cmake \
    ca-certificates \
    && rm -rf /var/lib/apt/lists/*

WORKDIR /src
COPY OpenWAM /src

RUN cmake -S . -B build -DCMAKE_BUILD_TYPE=Release \
    && cmake --build build --target OpenWAM -j"$(nproc)"

WORKDIR /work
ENTRYPOINT ["/src/build/bin/release/OpenWAM"]
