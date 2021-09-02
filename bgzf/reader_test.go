package bgzf

import (
	"bytes"
	"crypto/md5"
	"io"
	"io/ioutil"
	"os"
	"testing"
)

var block1hash = [16]byte{3, 56, 205, 200, 44, 174, 95, 168, 153, 248, 173, 154, 105, 114, 7, 32}
var block2hash = [16]byte{218, 54, 146, 68, 164, 54, 39, 51, 148, 255, 120, 50, 109, 214, 29, 2}

func TestRead(t *testing.T) {
	b := NewBlock()
	r := NewBlockReader("testdata/test.bam")
	var err error

	b1Offset := getHeaderOffset(r) + 1
	err = r.ReadBlock(b)
	if err != nil || md5.Sum(b.Bytes()) != block1hash {
		t.Error("problem reading BGZF file")
	}

	b2Offset := getHeaderOffset(r) + 1 + b1Offset
	err = r.ReadBlock(b)
	if err != nil || md5.Sum(b.Bytes()) != block2hash {
		t.Error("problem reading BGZF file")
	}

	err = r.ReadBlock(b)
	if err != io.EOF {
		t.Error("problem reading BGZF file")
	}

	err = r.Close()
	if err != nil {
		t.Error("problem reading BGZF file")
	}

	// Seek tests
	b = NewBlock()
	r = NewBlockReader("testdata/test.bam")

	_, err = r.Seek(b1Offset, io.SeekStart)
	if err != nil {
		t.Error("problem with BGZF seek")
	}
	err = r.ReadBlock(b)
	if err != nil || md5.Sum(b.Bytes()) != block2hash {
		t.Error("problem with BGZF seek")
	}

	_, err = r.Seek(b2Offset, io.SeekStart)
	if err != nil {
		t.Error("problem with BGZF seek")
	}
	err = r.ReadBlock(b)
	if err != io.EOF {
		t.Error("problem with BGZF seek")
	}

	_, err = r.Seek(0, io.SeekStart)
	if err != nil {
		t.Errorf("problem with BGZF seek")
	}
	err = r.ReadBlock(b)
	if err != nil || md5.Sum(b.Bytes()) != block1hash {
		t.Error("problem with BGZF seek")
	}

	err = r.Close()
	if err != nil {
		t.Error("problem reading BGZF file")
	}
}

func TestByteReadWrite(t *testing.T) {
	actual, _ := os.Open("testdata/test.sam")
	var actualBytes, readBytes bytes.Buffer
	io.Copy(&actualBytes, actual)
	actual.Close()

	// write actual bytes to bgzf writer
	tmpFile, _ := ioutil.TempFile("", "")
	writer := NewWriter(tmpFile)
	writer.Write(actualBytes.Bytes())
	err := writer.Close()
	if err != nil {
		t.Errorf("problem with bgzf byte writer")
	}
	tmpFile.Close()

	// read written bytes
	testFile := NewReader(tmpFile.Name())
	io.Copy(&readBytes, testFile)
	err = testFile.Close()
	if err != nil {
		t.Errorf("problem with bgzf byte reader")
	}

	if !bytes.Equal(actualBytes.Bytes(), readBytes.Bytes()) {
		t.Errorf("problem with bgzf byte reader")
	}

	os.Remove(tmpFile.Name())
}
