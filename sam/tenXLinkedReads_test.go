package sam

import (
	"testing"
)

func Test10xBarcode(t *testing.T) {
	TenXPrettyPrint("testdata/RABS.10x.sam")
}

func TestWriting10xBam(t *testing.T) {
	SwitchBxTagChannel("testdata/RABS.10x.sam", "testdata/RABS.10xProcessed.sam", 8)
}
