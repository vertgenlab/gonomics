package sam

import (
	"fmt"
	"testing"
)

func TestTags(t *testing.T) {
	var err error
	data, _ := Read("testdata/unmapped.bam")
	for i := range data {
		err = ParseExtra(&data[i])
		if err != nil {
			t.Error("problem parsing tags in ParseExtra")
		}
	}

	expected, _ := Read("testdata/unmapped.sam")

	for i := range data {
		if data[i].Extra != expected[i].Extra {
			fmt.Println(data[i].Extra)
			fmt.Println(expected[i].Extra)
			t.Error("problem parsing tags in ParseExtra")
		}
	}
}
