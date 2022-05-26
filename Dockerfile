FROM golang:1.18-buster as builder

ENV TZ=America/New_York
RUN ln -snf /usr/share/zoneinfo/$TZ /etc/localtime && echo $TZ > /etc/timezone

# Update linux image which contains golang-1.18
RUN set -x && apt-get update && DEBIAN_FRONTEND=noninteractive apt-get install -y \
  ca-certificates && \
  rm -rf /var/lib/apt/lists/*

# To build docker container image: run the following command:
#
### docker build . -t gonomics
# 
# Use gonomics programs with:
#
### docker run gonomics gsw -h
#
# All binary executables should exist in /bin if built sucessfully, to check, run:
### docker run -it gonomics ls /go/bin
#

# Set up working directory and add all golang files to container
WORKDIR /src/github.com/vertgenlab/gonomics
ADD . ./

# Download all gonomics dependencies and modules
RUN go mod download \
  && go mod vendor \
  && go mod verify

# Install all binary executables into $GOBIN
RUN go install ./...

