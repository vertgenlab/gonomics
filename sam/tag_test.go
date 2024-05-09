package sam

import (
	"fmt"
	"github.com/vertgenlab/gonomics/exception"
	"io"
	"os"
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

func TestBamWriteTags(t *testing.T) {
	file, err := os.CreateTemp("", "")
	data, header := Read("testdata/unmapped.sam")
	w := NewBamWriter(file, header)
	for i := range data {
		WriteToBamFileHandle(w, data[i], 0)
	}
	err = w.Close()
	exception.PanicOnErr(err)
	file.Close()

	var afterWrite []Sam
	br, _ := OpenBam(file.Name())
	for err != io.EOF {
		var curr Sam
		_, err = DecodeBam(br, &curr)
		ParseExtra(&curr)
		if err == nil {
			afterWrite = append(afterWrite, curr)
		}
	}
	err = br.Close()
	exception.PanicOnErr(err)

	for i := range data {
		if data[i].Extra != afterWrite[i].Extra {
			fmt.Println(data[i].Extra)
			fmt.Println(afterWrite[i].Extra)
			t.Error("problem writing tags in bam files")
		}
	}
}
