package fileio

import (
	"bufio"
	"compress/gzip"
	"fmt"
	"net/http"
	"strings"

	"github.com/vertgenlab/gonomics/exception"
)

// EasyHttp will fetch data from files uploaded to an internet server and stream data into EasyReader.
func EasyHttp(url string) *EasyReader {
	answer := EasyReader{}
	resp, err := http.Get(url)
	exception.PanicOnErr(err)
	if strings.HasSuffix(url, ".gz") {
		answer.internalGzip, err = gzip.NewReader(resp.Body)
		exception.PanicOnErr(err)
		answer.BuffReader = bufio.NewReader(answer.internalGzip)
	} else {
		answer.BuffReader = bufio.NewReader(resp.Body)
		answer.internalGzip = nil
	}
	return &answer
}

// CatUrl will process a url link and print it out to stdout.
func CatUrl(url string) string {
	var answer string
	reader := EasyOpen(url)
	for i, done := EasyNextLine(reader); !done; i, done = EasyNextLine(reader) {
		answer += fmt.Sprintf("%s\n", i)
	}
	return answer
}
