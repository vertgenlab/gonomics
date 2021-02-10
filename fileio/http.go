package fileio

import (
	"bufio"
	"compress/gzip"
	"fmt"
	"net/http"
	"strings"
)

// EasyHttp will fetch data from files uploaded to an internet server and stream data into EasyReader
func EasyHttp(url string) *EasyReader {
	answer := EasyReader{}
	resp, err := http.Get(url)
	panicOnErr(err)
	if strings.HasSuffix(url, ".gz") {
		answer.internalGzip, err = gzip.NewReader(resp.Body)
		panicOnErr(err)
		answer.BuffReader = bufio.NewReader(answer.internalGzip)
	} else {
		answer.BuffReader = bufio.NewReader(resp.Body)
		answer.internalGzip = nil
	}
	return &answer
}

// CatUrl will process a url link and print it out to stdout.
func CatUrl(url string) {
	reader := EasyOpen(url)
	for i, done := EasyNextLine(reader); !done; i, done = EasyNextLine(reader) {
		fmt.Printf("%s\n", i)
	}
}
