package fileio

import (
	"bufio"
	"bytes"
	"compress/gzip"
	"fmt"
	"github.com/vertgenlab/gonomics/common"
	"net/http"
	"strings"
)

// HttpReader will fetch data from files uploaded to an internet server and stream data into an io.Reader interface.
func HttpReader(url string) *SimpleReader {
	resp, err := http.Get(url)
	common.ExitIfError(err)
	var answer SimpleReader = SimpleReader{
		Buffer: &bytes.Buffer{},
		line:   make([]byte, defaultBufSize),
	}
	if strings.HasSuffix(url, ".gz") {
		gzipReader, err := gzip.NewReader(resp.Body)
		common.ExitIfError(err)
		answer.Reader = bufio.NewReader(gzipReader)
	} else {
		answer.Reader = bufio.NewReader(resp.Body)
	}
	return &answer
}

// VimUrl is a basic function to procress a url link and print it out to stdout.
func VimUrl(url string) {
	reader := HttpReader(url)
	for i, err := ReadLine(reader); !err; i, err = ReadLine(reader) {
		fmt.Printf("%s\n", i.String())
	}
}

// ReadUrl is a basic function to fetch data and read data from internet files and return a slice of strings for each line.
func ReadUrl(url string) []string {
	var answer []string
	reader := HttpReader(url)
	for i, err := ReadLine(reader); !err; i, err = ReadLine(reader) {
		answer = append(answer, i.String())
	}
	return answer
}
