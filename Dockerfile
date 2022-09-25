FROM debian:bullseye-slim as downloader

RUN apt-get update && \ 
    apt-get install -y build-essential \
                       cmake \
                       git \
                       ssh \
                       perl \
                       python3 \
                       autoconf \
                       automake \
                       libtool \
                       uuid-dev \
                       gcovr

# Download build and install Catch2 framework
WORKDIR /root
RUN git clone https://github.com/catchorg/Catch2 && \
    cd Catch2 && \
    git checkout v3.1.0 && \
    cmake -B build -S . -DCMAKE_INSTALL_PREFIX=/usr/local && \
    cmake --build build/ --target install && \
    ldconfig

# Copy all from downloader image to strip unnecessary files from final image
FROM debian:bullseye-slim

COPY --from=downloader /bin/ /bin
COPY --from=downloader /etc/ /etc
COPY --from=downloader /lib/ /lib
COPY --from=downloader /lib64/ /lib
COPY --from=downloader /opt/ /opt
COPY --from=downloader /sys/ /sys
COPY --from=downloader /usr/ /usr

RUN ldconfig

RUN adduser --quiet user

USER user

CMD cmake -B ./build -S . -DCMAKE_BUILD_TYPE=Release \
                          -DCMAKE_INSTALL_PREFIX=/usr/local \
                          -DENABLE_TESTS=ON && \
    cmake --build ./build -j`nproc`