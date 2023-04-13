all: install test lint

build:
	docker build . -t gonomics

clean:
	go fmt ./...
	golangci-lint run ./... --fix

install:
	go mod download && go mod verify
	go install ./...
	go install github.com/golangci/golangci-lint/cmd/golangci-lint@v1.47.3

lint: test
	golangci-lint run ./...

list:
	@grep '^[^#[:space:]].*:' Makefile

test:
	go test ./...
