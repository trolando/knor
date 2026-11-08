# knor for SYNTCOMP's parity-game realizability track, as the competition runs
# it since 2025 (Docker images on the TGCC cluster, see README.md): Debian
# bookworm with runlim, and in /home/solvers one script per configuration that
# takes the eHOA file as its only argument.
#
#     docker build --platform linux/amd64 -t knor-syntcomp .
#     docker run --rm -v "$PWD/examples:/in:ro" -v "$PWD/test/expected_verdicts.txt:/expected:ro" \
#         knor-syntcomp /usr/local/bin/syntcomp-test /in /expected

FROM gcc:13.3.0-bookworm AS build

RUN apt-get update \
    && apt-get install -y --no-install-recommends cmake libboost-filesystem-dev libboost-iostreams-dev \
        libboost-random-dev libboost-regex-dev libboost-system-dev \
    && rm -rf /var/lib/apt/lists/*

# runlim 2.0.0rc12, the version in the competition's images
RUN git clone https://github.com/arminbiere/runlim.git /runlim \
    && git -C /runlim checkout 3076e4406466ee14c3b6decc86c31b25c570cdc6 \
    && cd /runlim && ./configure.sh && make runlim

# The binary must run on the competition's machines, not only on the one that
# builds it: no -march=native (Oink has an option for it, Lace and Sylvan add
# it unconditionally).
COPY . /knor
RUN sed -i 's/-march=native//' /knor/libs/lace/CMakeLists.txt /knor/libs/sylvan/src/CMakeLists.txt \
    && cmake -S /knor -B /build -DCMAKE_BUILD_TYPE=Release -DOINK_BUILD_NATIVE=OFF \
    && cmake --build /build --target knor -j"$(nproc)"

FROM debian:bookworm

COPY --from=build /runlim/runlim /usr/bin/runlim
COPY --from=build /build/knor /usr/local/bin/knor
COPY syntcomp/solvers/ /home/solvers/
COPY syntcomp/test.sh /usr/local/bin/syntcomp-test

ENTRYPOINT ["/bin/bash"]
