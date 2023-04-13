install:
	go mod download && go mod verify
	go install ./...
	go install github.com/golangci/golangci-lint/cmd/golangci-lint@v1.47.3

test:
	go test ./...

lint:
	golangci-lint run ./...

clean:
	golangci-lint run ./... --fix

build:
	docker build . -t gonomics
