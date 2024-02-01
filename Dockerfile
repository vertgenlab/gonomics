FROM golang:1.20.0-alpine3.16

# Build docker image with the following command:
# docker build . -t gonomics

# Example: Use gonomics programs with: docker run gonomics gsw -h

# Set up working directory and add all golang files to container
WORKDIR /src/github.com/vertgenlab/gonomics
ENV GOBIN=/bin

COPY go.mod go.sum ./
# Download all gonomics dependencies and modules
RUN go mod download && go mod verify

# Add gonomics repository files into container
COPY . .

# Install all binary executables into $GOBIN
RUN CGO_ENABLED=0 go install -v -ldflags='-w -s' ./...
