package sam

import (
	"fmt"
	"testing"
)

const bigFile string = "/Users/danielsnellings/Desktop/Ex_01_Results/final/Ex_01_01.grouped.bam.bai"

func TestReadBai(t *testing.T) {
	bai := ReadBai(bigFile)

	for _, bin := range bai.refs[50].bins {
		fmt.Sprintln(bin.id, bin.refStart, bin.refEnd)
	}
}
