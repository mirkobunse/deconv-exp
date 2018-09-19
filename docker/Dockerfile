FROM openjdk:8

MAINTAINER Mirko Bunse <mirko.bunse@cs.tu-dortmund.de>

# configuration
ARG group
ARG gid
ARG user
ARG uid

# create user
RUN mkdir -p /opt/home && \
    groupadd --gid $gid $group && \
    useradd  --gid $gid --uid $uid $user

# auxiliary tools
RUN apt-get update && apt-get upgrade -y && apt-get install -y \
    nano \
    less \
    hdf5-tools \
    imagemagick		# imagemagick also required by julia packages

# julia 0.6, following https://hub.docker.com/r/julialang/julia/
RUN mkdir -p /opt/julia-0.6.4 && \
    curl -s -L https://julialang-s3.julialang.org/bin/linux/x64/0.6/julia-0.6.4-linux-x86_64.tar.gz | tar -C /opt/julia-0.6.4 -x -z --strip-components=1 -f - && \
    ln -fs /opt/julia-0.6.4     /opt/julia && \
    ln -fs /opt/julia/bin/julia /bin/julia && \
    rm /opt/julia/etc/julia/juliarc.jl

# install whatever is needed by julia packages
RUN apt-get install -y \
    build-essential \
    make \
    libgtk-3-0 \
    gettext \
    libav-tools \
    && apt-get clean

# add startup routine
ADD juliarc.jl /opt/julia/etc/julia
ADD entrypoint.sh /opt
RUN chown --recursive $uid /opt

USER $uid
WORKDIR /home/$user/mt-exp

ENTRYPOINT ["/opt/entrypoint.sh"]
CMD ["bash"]

