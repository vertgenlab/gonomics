package bgzf

import (
	"crypto/md5"
	"io"
	"testing"
)

var block1hash = [16]byte{3, 56, 205, 200, 44, 174, 95, 168, 153, 248, 173, 154, 105, 114, 7, 32}
var block2hash = [16]byte{218, 54, 146, 68, 164, 54, 39, 51, 148, 255, 120, 50, 109, 214, 29, 2}

func TestRead(t *testing.T) {
	b := NewBlock()
	r := NewReader("testdata/test.bam")
	var err error

	err = r.ReadBlock(b)
	if err != nil || md5.Sum(b.Bytes()) != block1hash {
		t.Error("problem reading BGZF file")
	}

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
}

func TestSeek(t *testing.T) {
	b := NewBlock()
	r := NewReader("testdata/test.bam")
	var err error

	_, err = r.Seek(2299, io.SeekStart)
	if err != nil {
		t.Error("problem with BGZF seek")
	}
	err = r.ReadBlock(b)
	if err != nil || md5.Sum(b.Bytes()) != block2hash {
		t.Error("problem with BGZF seek")
	}

	_, err = r.Seek(2299+1217, io.SeekStart)
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
