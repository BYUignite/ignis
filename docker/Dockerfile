# Build the image: 
#     docker build -t ignis_img .
#
# Create and run in a container: 
#     docker run --name ignis_container -v /path/to/ignis/docker:/home --hostname docker -it ignis_img
# 
#     This is setup for an interactive shell. 
#     The folder ~/docker will be created in you home (host) if not already there and connected to /home in the docker container
#     This makes it easy to access and share files back and forth
#     add the --rm flag to delete the container on exit, like docker run --rm --name ...
#
# Connect to and existing container: 
#     docker start -ai ignis_container
# 
#     -a attaches, and -i is interactive
#
# When entering the docker container you will start in /app where ignis is installed.
# You can do "cp -r ignis /home" to copy the install to /home/ignis.
# This location is mounted to the host at /path/to/ignis/docker (in the docker run command above).
# You will see that /path/to/ignis/docker/ignis is now present on the host 
# and matches /home/ignis in the container.
# Inside the container you can run the ignis code, 
# and then view/plot/edit the files on your local host.
#
#################################################################################

FROM alpine:3.19.0 as buildit

RUN apk add --no-cache wget tar cmake make g++ lapack-dev vim

############## Build Sundials

WORKDIR /app/

RUN wget https://github.com/LLNL/sundials/releases/download/v6.7.0/sundials-6.7.0.tar.gz \
&&  tar -xvzf sundials-6.7.0.tar.gz \
&&  rm sundials-6.7.0.tar.gz \
&&  cd sundials-6.7.0 \
&&  mkdir build \
&&  cd build \
&&  cmake .. \
&&  make -j8 \
&&  make install \
&&  cd ../../ \
&&  rm -r sundials-6.7.0

############## Build Yaml

RUN apk add --no-cache git

RUN git clone https://github.com/jbeder/yaml-cpp.git -b 0.8.0 --depth=1 \
&&  cd yaml-cpp \
&&  mkdir build \
&&  cd build \
&&  cmake .. \
&&  make -j8 \
&&  make install \
&&  cd ../../ \
&&  rm -r yaml-cpp

############## Build Cantera

RUN apk add --no-cache scons python3 boost-dev hdf5-dev
RUN apk add --no-cache py3-setuptools

RUN git clone https://github.com/Cantera/cantera.git -b 3.0 --depth=1 \
&&  cd cantera \
&&  scons build -j8 python_package=n \
&&  scons install \
&&  cd ../ \
&&  rm -r cantera

############## Build radlib

RUN git clone https://github.com/BYUignite/radlib.git -b v1.3 --depth=1 \
&&  cd radlib \
&&  mkdir build \
&&  cd build \
&&  cmake --install-prefix /opt/radlib .. \
&&  make -j \
&&  make install \
&&  cd ../../ \
&&  rm -r radlib \
&&  export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/opt/radlib/lib

############## Build sootlib

RUN git clone https://github.com/BYUignite/sootlib.git -b v1.01 --depth=1 \
&&  cd sootlib \
&&  mkdir build \
&&  cd build \
&&  cmake --install-prefix /opt/sootlib .. \
&&  make -j \
&&  make install \
&&  cd ../../ \
&&  rm -r sootlib \
&&  export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/opt/sootlib/lib

############## Build ignis

RUN git clone https://github.com/BYUignite/ignis.git --depth=1 \
&&  cd ignis \
&&  mkdir build \
&&  cd build \
&&  cmake .. \
&&  make -j \
&&  make install \
&&  cd ../

############## Set environment variables for containers that use this image

ENV LD_LIBRARY_PATH=/opt/sootlib/lib:/opt/radlib/lib
