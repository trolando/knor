# Knor for SYNTCOMP. Ship sources and build dependencies so scripts/compile.sh
# can compile on the machine that will run the solver.
#
#     docker build -t knor .
#     docker run --rm --entrypoint /bin/bash -v knor-build:/opt/knor/build \
#         knor ./scripts/compile.sh
#     docker run --rm -i -v knor-build:/opt/knor/build:ro \
#         knor --paritize --auto-encode --compress -a < spec.ehoa

FROM debian:trixie-slim

ARG DEBIAN_FRONTEND=noninteractive

RUN apt-get update \
    && apt-get install -y --no-install-recommends \
        g++ cmake make libboost-random-dev ca-certificates \
    && rm -rf /var/lib/apt/lists/*

WORKDIR /opt/knor

# Initialize submodules before building the image; compilation needs their
# complete source trees and does not fetch anything at runtime.
COPY . .
RUN chmod +x scripts/compile.sh scripts/knor.sh scripts/knor-synt.sh

ENTRYPOINT ["/opt/knor/scripts/knor.sh"]
