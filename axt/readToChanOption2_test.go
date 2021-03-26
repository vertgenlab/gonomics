package axt

import (
	"fmt"
	"testing"
)

func TestGoReadToChanWithHeader(t *testing.T) {
	dataChan, header := GoReadToChanWithHeader("testdata/chrM_gasacu1.axt")
	fmt.Sprintln(header)
	for record := range dataChan {
		fmt.Sprintln(record)
	}
}
