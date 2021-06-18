package sam

import (
	"fmt"
	"testing"
)

const bamTestfile string = "../bgzf/testdata/test.bam"

func TestReadBam(t *testing.T) {
	_, header, chroms := OpenBam(bamTestfile)

	if len(header.Chroms) != len(chroms) {
		t.Errorf("plain text and bam header chroms do not match")
	}

	for i := range chroms {
		if chroms[i] != header.Chroms[i] {
			fmt.Println(len(chroms[i].Name), len(header.Chroms[i].Name))
			t.Errorf("plain text and bam header chroms do not match")
		}
	}
}
