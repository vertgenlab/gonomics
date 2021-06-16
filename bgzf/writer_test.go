package bgzf

import (
	"bytes"
	"compress/gzip"
	"encoding/hex"
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
		_, WriteErr := w.Write(b.Bytes())
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

	//r = NewReader("testdata/tmp.bam")
	//fmt.Println(r.decompressor.Header.Extra)
	//err = r.ReadBlock(b)
	//fmt.Println(b.Len())
	////fmt.Println(string(b.Bytes()))
	//
	//fmt.Println(r.decompressor.Header.Extra)
	//err = r.ReadBlock(b)
	//fmt.Println(b.Len())
	////fmt.Println(string(b.Bytes()))
	//
	//fmt.Println(r.decompressor.Header.Extra)
	//err = r.ReadBlock(b)
	//fmt.Println(b.Len())
	////fmt.Println(string(b.Bytes()))

	if !equalFiles("testdata/test.bam", "testdata/tmp.bam") {
		t.Error("problem writing BGZF file")
	}

	//err = os.Remove("testdata/tmp.bam")
	//if err != nil {
	//	log.Panic(err)
	//}
}

func equalFiles(aname, bname string) bool {
	ablk := NewBlock()
	bblk := NewBlock()
	a := NewReader(aname)
	b := NewReader(bname)
	var adata, bdata []byte
	var aErr, bErr error
	for {
		fmt.Println("A: ", a.decompressor.Header.Extra)
		fmt.Println("B: ", b.decompressor.Header.Extra)
		aErr = a.ReadBlock(ablk)
		bErr = b.ReadBlock(bblk)
		fmt.Println()

		fmt.Println(aErr, bErr)
		if aErr != bErr {
			return false
		}

		if aErr == bErr && aErr == io.EOF {
			return true
		}

		fmt.Println(ablk.Len(), bblk.Len())
		if ablk.Len() != bblk.Len() {
			return false
		}
		adata = ablk.Bytes()
		bdata = ablk.Bytes()
		for i := range adata {
			if adata[i] != bdata[i] {
				return false
			}
		}
	}
}

func TestCompress(t *testing.T) {
	var a bytes.Buffer
	zr, _ := gzip.NewWriterLevel(&a, gzip.DefaultCompression)
	zr.Reset(&a)
	fmt.Println(zr.Header)
	zr.Header.Extra = []byte{66,67,2,0,27,0}
	zr.Write([]byte{})
	zr.Close()
	fmt.Println("empty compress is size: ", a.Len())
	fmt.Println(hex.EncodeToString(a.Bytes()))
}
