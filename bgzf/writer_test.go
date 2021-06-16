package bgzf

import (
	"bytes"
	"compress/gzip"
	"crypto/md5"
	"encoding/binary"
	"fmt"
	"io"
	"log"
	"os"
	"testing"
)

func TestWrite(t *testing.T) {
	b := NewBlock()
	r := NewReader("testdata/test.bam")
	tmpfile, err := os.Create("testdata/tmp.bam")
	if err != nil {
		log.Panic(err)
	}
	w := NewWriter(tmpfile)

	for err = r.ReadBlock(b); err != io.EOF; err = r.ReadBlock(b) {
		n, WriteErr := w.Write(b.Bytes())
		if n != b.Len() {
			t.Error("write size differs")
		}
		if WriteErr != nil {
			t.Error("problem writing BGZF file")
		}
	}

	err = r.Close()
	if err != nil {
		log.Panic(err)
	}
	err = w.Close()
	if err != nil {
		t.Error("problem writing BGZF file: ", err)
	}
	err = tmpfile.Close()
	if err != nil {
		log.Panic(err)
	}

	r = NewReader("testdata/tmp.bam")

	if len(r.decompressor.Header.Extra) < 6 {
		t.Error("bgzf header was not written")
	}
	b1Offset := getHeaderOffset(r) + 1
	err = r.ReadBlock(b)
	if err != nil || md5.Sum(b.Bytes()) != block1hash {
		t.Error("problem reading BGZF file")
	}

	if len(r.decompressor.Header.Extra) < 6 {
		t.Error("bgzf header was not written")
	}
	b2Offset := getHeaderOffset(r) + 1 + b1Offset
	err = r.ReadBlock(b)
	if err != nil || md5.Sum(b.Bytes()) != block2hash {
		t.Error("problem reading BGZF file")
	}

	if len(r.decompressor.Header.Extra) < 6 {
		t.Error("bgzf header was not written")
	}
	err = r.ReadBlock(b)
	if err != io.EOF {
		t.Error("problem reading BGZF file")
	}

	//Seek tests
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

	tmpfile.Close()
	os.Remove("testdata/tmp.bam")
}

func TestSimpleCompress(t *testing.T) {
	var a bytes.Buffer
	var data []byte = []byte("hello\n")
	r := NewWriter(&a)
	n, err := r.Write(data)
	if err != nil {
		fmt.Printf("input %d byte, wrote %d bytes\n", len(data), n)
		t.Error("trouble writing: ", err)
	}
	err = r.Close()
	if err != nil {
		t.Error("trouble writing: ", err)
	}

	nr, err := gzip.NewReader(&a)
	if err != nil {
		t.Error(err)
	}

	var b bytes.Buffer
	_, err = io.Copy(&b, nr)
	if err != nil {
		t.Error(err)
	}
	if b.String() != "hello\n" {
		t.Error("problem with bgzf compression")
	}
}

func getHeaderOffset(r Reader) int64 {
	extra := r.decompressor.Header.Extra
	return int64(binary.LittleEndian.Uint16(extra[4:6]))
}
