all: install test lint

build:
	docker build . -t gonomics

clean:
	go fmt ./...
	gofmt -s -w ./*/*.go
	gofmt -s -w ./*/*/*.go
	gofmt -s -w ./cmd/*/*.go

install:
	go mod download && go mod verify
	go install ./...

list:
	@grep '^[^#[:space:]].*:' Makefile

test:
	go test -covermode=atomic ./...
