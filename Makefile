all: install test lint

build:
	docker build . -t gonomics

clean:
	go fmt ./...
	gofmt -s -w ./*/*.go
	gofmt -s -w ./*/*/*.go
	gofmt -s -w ./cmd/*/*.go
	golangci-lint run ./... --fix

install:
	go mod download && go mod verify
	go install ./...
	curl -sSfL https://raw.githubusercontent.com/golangci/golangci-lint/master/install.sh | sh -s -- -b $(go env GOPATH)/bin v1.52.2

lint: test
	golangci-lint run ./...

list:
	@grep '^[^#[:space:]].*:' Makefile

test:
	go test -covermode=atomic ./...
