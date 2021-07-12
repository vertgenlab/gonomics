package sam

import (
	"bytes"
	"fmt"
	"testing"
)

func TestBamWriter(t *testing.T) {
	var err error
	actual := new(bytes.Buffer)
	bw := NewBamWriter(actual, Header{})
	err = bw.Close()
	if err != nil {
		t.Error(err)
	}

	fmt.Sprintln(bw)
	fmt.Println(actual)
}
