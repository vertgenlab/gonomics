package sam

import (
	"fmt"
	"github.com/vertgenlab/gonomics/fileio"
	"testing"
	"unsafe"
)

func TestReadHeaderBytes(t *testing.T) {
	br := fileio.NewByteReader("testdata/test.sam")
	header := ReadHeaderBR(br)

	printMetadata(header.Metadata)
	var a Aln
	fmt.Println("Sam Size is: ", unsafe.Sizeof(a))
}

func printMetadata(m Metadata) {
	fmt.Printf("Version: %s\n", m.Version)
	fmt.Printf("Sort Order: %s\n", m.SortOrder)
	fmt.Printf("Grouping: %s\n", m.Grouping)

	for key, val := range m.AllTags {
		fmt.Printf("Tag: %s\n", string(key[:]))
		var taglineString string
		for _, tagline := range val {
			var substring string
			for subkey, subval := range tagline {
				substring += fmt.Sprintf("\t%s:%s", string(subkey[:]), subval)
			}
			taglineString += substring + "\n"
		}
		fmt.Println(taglineString)
	}
}
