FROM golang:1.20.3-alpine

# To build docker container image: run the following command: docker build . -t gonomics
# Use gonomics programs with: docker run gonomics gsw -h
# All binary executables should exist in /bin if built sucessfully, to check, run:
# docker run -it gonomics ls /go/bin

# Set up working directory and add all golang files to container
WORKDIR /src/github.com/vertgenlab/gonomics
COPY go.mod go.sum ./
# Download all gonomics dependencies and modules
RUN go mod download && go mod verify
# Add gonomics repository files into container
COPY . .
# Install all binary executables into $GOBIN
RUN go install ./...
