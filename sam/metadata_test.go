package sam

import (
	"fmt"
	"github.com/vertgenlab/gonomics/chromInfo"
	"github.com/vertgenlab/gonomics/fileio"
	"testing"
)

func TestReadHeaderBytes(t *testing.T) {
	br := fileio.EasyOpen("testdata/small.sam")
	header := ReadHeader(br)

	expectedChrom := chromInfo.ChromInfo{Name: "ref", Size: 45, Order: 0}
	if len(header.Chroms) != 1 || header.Chroms[0] != expectedChrom {
		t.Error("problem with reading sam header")
	}

	if header.Metadata.Version != "1.6" {
		t.Error("problem with reading sam header")
	}

	expectedComment := "@CO\tthis is a comment line and\tmay\thave\ttabs"
	if len(header.Metadata.Comments) != 1 || header.Metadata.Comments[0] != expectedComment {
		t.Error("problem with reading sam header")
	}

	if header.Metadata.Grouping != Reference {
		t.Error("problem with reading sam header")
	}

	if len(header.Metadata.SortOrder) != 1 ||
		header.Metadata.SortOrder[0] != Coordinate {
		t.Error("problem with reading sam header")
	}

	tagmap := header.Metadata.AllTags
	if tagmap[[2]byte{'T', 'S'}][0][[2]byte{'S', 'S'}] != "this is a test" {
		t.Error("problem with reading sam header")
	}

	if tagmap[[2]byte{'T', 'S'}][0][[2]byte{'X', 'Q'}] != "777" {
		t.Error("problem with reading sam header")
	}
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
