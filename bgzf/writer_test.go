package bgzf

import (
	"bytes"
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
		_, WriteErr := w.Write(b.Bytes())
		if WriteErr != nil {
			t.Error("problem writing BGZF file")
		}
	}

	err = r.Close()
	if err != nil {
		log.Panic(err)
	}
	err = tmpfile.Close()
	if err != nil {
		log.Panic(err)
	}
	err = w.Close()
	if err != nil {
		t.Error("problem writing BGZF file: ", err)
	}

	if !equalFiles("testdata/test.bam", "testdata/tmp.bam") {
		t.Error("problem writing BGZF file")
	}

	err = os.Remove("testdata/tmp.bam")
	if err != nil {
		log.Panic(err)
	}
}

func equalFiles(aname, bname string) bool {
	a, _ := os.Open(aname)
	b, _ := os.Open(bname)

	var abuf, bbuf bytes.Buffer
	aRead, aErr := abuf.ReadFrom(a)
	bRead, bErr := bbuf.ReadFrom(b)

	if aRead != bRead || aErr != bErr || aErr != nil {
		return false
	}

	var aval, bval byte
	for abuf.Len() > 0 {
		aval, _ = abuf.ReadByte()
		bval, _ = bbuf.ReadByte()
		if aval != bval {
			return false
		}
	}
	return true
}
