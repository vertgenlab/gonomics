package sam

import (
	"bytes"
	"fmt"
	"io"
	"io/ioutil"
	"os"
	"strings"
	"testing"
)

func TestBamWriter(t *testing.T) {
	var err error
	data, header := Read("testdata/small.sam")

	actual := new(bytes.Buffer)
	bw := NewBamWriter(actual, header)

	for i := range data {
		WriteToBamFileHandle(bw, data[i], 0)
	}

	err = bw.Close()
	if err != nil {
		t.Error(err)
	}

	f, err := ioutil.TempFile("testdata", "tmp*.bam")
	if err != nil {
		t.Error(err)
	}
	io.Copy(f, actual)
	f.Close()

	actualData, actualHeader := GoReadToChan(f.Name())

	if strings.Join(actualHeader.Text, "\n") != strings.Join(header.Text, "\n") {
		t.Error("issue with writing header")
	}

	var i int
	for actualRecord := range actualData {
		if actualRecord.String() != data[i].String() {
			fmt.Println(actualRecord)
			fmt.Println(data[i])
			fmt.Println()
			t.Error("issue with writing sam as bam")
		}
		i++
	}

	os.Remove(f.Name())
}
