package sam

import (
	"fmt"
	"io"
	"testing"
)

const bamTestfile string = "../bgzf/testdata/test.bam"

func TestReadBam(t *testing.T) {
	r, header := OpenBam(bamTestfile)

	if len(header.Chroms) != len(r.refs) {
		t.Errorf("plain text and bam header chroms do not match")
	}

	for i := range r.refs {
		if r.refs[i] != header.Chroms[i] {
			t.Errorf("plain text and bam header chroms do not match")
		}
	}

	var err error
	for err != io.EOF {
		_, _, err = NextBam(r)
		fmt.Println("read bam")
	}

}
