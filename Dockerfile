FROM ubuntu:22.04 AS builder
RUN apt-get update && apt-get install -y build-essential cmake git
WORKDIR /build
COPY scpsol/ scpsol/
COPY data/ data/
WORKDIR /build/scpsol
RUN cmake -B build -DCMAKE_BUILD_TYPE=Release \
    && cmake --build build -j$(nproc)

FROM ubuntu:22.04
COPY --from=builder /build/scpsol/build/scpsol /usr/local/bin/scpsol
COPY data/ /data/
ENTRYPOINT ["scpsol"]
CMD ["/data/scp_demo00.txt"]
