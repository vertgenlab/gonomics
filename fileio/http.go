package fileio

import (
	"bufio"
	"compress/gzip"
	"fmt"
	"github.com/vertgenlab/gonomics/common"
	"net/http"
	"strings"
)

// EasyHttp will fetch data from files uploaded to an internet server and stream data into EasyReader
func EasyHttp(url string) *EasyReader {
	answer := EasyReader{}
	resp, err := http.Get(url)
	common.ExitIfError(err)
	if strings.HasSuffix(url, ".gz") {
		answer.internalGzip, err = gzip.NewReader(resp.Body)
		common.ExitIfError(err)
		answer.BuffReader = bufio.NewReader(answer.internalGzip)
	} else {
		answer.BuffReader = bufio.NewReader(resp.Body)
		answer.internalGzip = nil
	}
	return &answer
}

// CatUrl is a basic function to procress a url link and print it out to stdout.
func CatUrl(url string) {
	reader := EasyOpen(url)
	for i, done := EasyNextLine(reader); !done; i, done = EasyNextLine(reader) {
		fmt.Printf("%s\n", i)
	}
}
