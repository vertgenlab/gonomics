package sam

import (
	"testing"
)

const testfile string = "../bgzf/testdata/test.bam.bai"

func TestReadBai(t *testing.T) {
	bai := ReadBai(testfile)

	if len(bai.refs) != 84 {
		t.Errorf("problem reading bai file")
	}

	tRef := bai.refs[0]

	if tRef.head.children[0].id != 4681 {
		t.Errorf("problem reading bai file")
	}

	if tRef.bins[0].parent.id != tRef.head.id {
		t.Errorf("problem reading bai file")
	}

	// TODO more and better tests, especially of the tree structure
	// making these tests should be easier once a bam reader is
	// implemented so we can test seek behavior and virtual offsets.
}
