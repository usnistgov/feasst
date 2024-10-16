# Configuration Docker (Rocky 9)
# sudo dnf config-manager --add-repo https://download.docker.com/linux/centos/docker-ce.repo
# sudo dnf install docker-ce docker-ce-cli containerd.io
# sudo systemctl start docker
# sudo systemctl status docker
# sudo systemctl enable docker
# sudo docker build .

FROM ubuntu:24.04

RUN apt-get update

RUN apt-get install -y python3 python3-pip wget tar cmake g++

RUN  wget https://github.com/usnistgov/feasst/archive/refs/tags/v0.25.2.tar.gz \
  && tar -xf v0.25.2.tar.gz --transform s/-/-v/ \
  && mkdir feasst-v0.25.2/build; cd feasst-v0.25.2/build \
  && cmake .. \
  && make install -j$(nproc)

RUN ./feasst-v0.25.2/build/bin/fst < ./feasst-v0.25.2/tutorial/example.txt
