package sam

import (
	"bytes"
	"errors"
	"fmt"
	"io"
	"os"
	"strings"
	"testing"
)

func TestBamWriter(t *testing.T) {
	var err error = testFile("testdata/small.sam")
	if err != nil {
		t.Error(err)
	}
}

func TestBamWriterAuxTags(t *testing.T) {
	var err error = testFile("testdata/auxTagTest.sam")
	if err != nil {
		t.Error(err)
	}
}

func testFile(file string) error {
	var err error
	data, header := Read(file)

	actual := new(bytes.Buffer)
	bw := NewBamWriter(actual, header)

	for i := range data {
		WriteToBamFileHandle(bw, data[i], 0)
	}

	err = bw.Close()
	if err != nil {
		return err
	}

	f, err := os.CreateTemp("testdata", "tmp*.bam")
	if err != nil {
		return err
	}
	io.Copy(f, actual)
	f.Close()

	actualData, actualHeader := GoReadToChan(f.Name())

	if strings.Join(actualHeader.Text, "\n") != strings.Join(header.Text, "\n") {
		return errors.New("issue with writing header")
	}

	var i int
	for actualRecord := range actualData {
		if actualRecord.String() != data[i].String() {
			return fmt.Errorf("issue with writing sam as bam, record = %s, expected = %s ", actualRecord.String(), data[i].String())
		}
		i++
	}

	os.Remove(f.Name())
	return nil
}
