FROM ubuntu:18.04

RUN apt-get update

# miniconda3
RUN apt-get install wget -y && \
    wget http://repo.continuum.io/miniconda/Miniconda3-4.7.12.1-Linux-x86_64.sh -O /root/miniconda.sh && \
    chmod 0755 /root/miniconda.sh && \
    /root/miniconda.sh -b -p /root/miniconda3 && \
    rm /root/miniconda.sh && \
    apt-get remove wget -y
ENV PATH="/root/miniconda3/bin:$PATH"

# feasst
RUN apt-get install -y g++ cmake swig git && \
    git clone https://github.com/usnistgov/feasst.git && \
    mkdir feasst/build && cd feasst/build && \
    cmake -DUSE_SWIG=ON -DSET_PYTHON_PATH=ON -DPYTHON_INCLUDE_DIR=$HOME/miniconda3/include/python3.7m -DPYTHON_LIBRARY=$HOME/miniconda3/lib/libpython3.7m.so .. && \
    make _feasst -j$(nproc) && \
    make install -j$(nproc) && \
    apt-get remove -y cmake swig git && \
    apt-get autoremove -y && \
    pip install jupyter matplotlib scipy pandas
