package vcf

import (
	"github.com/vertgenlab/gonomics/chromInfo"
	"testing"
)

func expectedHeader() Header {
	var h Header
	h.FileFormat = "VCFv4.3"
	h.Chroms = map[string]chromInfo.ChromInfo{
		"chrA": {Name: "chrA", Size: 10, Order: 0},
		"chrB": {Name: "chrB", Size: 20, Order: 1}}
	h.Filter = map[string]FilterHeader{
		"FilterC": {Id: "FilterC", Description: "Sea"},
		"FilterD": {Id: "FilterD", Description: "Dea"},
	}
	h.Info = map[string]InfoHeader{
		"InfoA":      {Key: Key{Id: "InfoA", Number: "1", DataType: Integer, IsFormat: false}, Description: "Aye"},
		"InfoB":      {Key: Key{Id: "InfoB", Number: "R", DataType: Float, IsFormat: false}, Description: "Bee"},
		"InfoChar":   {Key: Key{Id: "InfoChar", Number: "2", DataType: Character, IsFormat: false}, Description: "Gee"},
		"InfoString": {Key: Key{Id: "InfoString", Number: "A", DataType: String, IsFormat: false}, Description: "Ach"},
		"InfoFlag":   {Key: Key{Id: "InfoFlag", Number: "0", DataType: Flag, IsFormat: false}, Description: "Eye"},
	}
	h.Format = map[string]FormatHeader{
		"GT":      {Key: Key{Id: "GT", Number: "1", DataType: Integer, IsFormat: true}, Description: "Eey"},
		"FormatF": {Key: Key{Id: "FormatF", Number: "1", DataType: Integer, IsFormat: true}, Description: "Eff"},
		"FormatJ": {Key: Key{Id: "FormatJ", Number: "2", DataType: Character, IsFormat: true}, Description: "Jay"},
		"FormatK": {Key: Key{Id: "FormatK", Number: ".", DataType: String, IsFormat: true}, Description: "Kay"},
	}
	h.Samples = map[string]int{
		"SampleA": 0,
		"SampleB": 1,
	}
	return h
}

func TestHeaderReading(t *testing.T) {
	_, aH := Read("testdata/headerTest.vcf") // actual header
	eH := expectedHeader()                   // expected header

	if aH.FileFormat != eH.FileFormat {
		t.Errorf("problem with header file format")
	}

	if len(aH.Chroms) != len(eH.Chroms) ||
		len(aH.Info) != len(eH.Info) ||
		len(aH.Filter) != len(eH.Filter) ||
		len(aH.Format) != len(eH.Format) ||
		len(aH.Samples) != len(eH.Samples) {
		t.Errorf("problem with header maps")
	}

	for key := range aH.Chroms {
		if aH.Chroms[key] != eH.Chroms[key] {
			t.Errorf("problem with header chroms")
		}
	}
	for key := range aH.Info {
		if aH.Info[key] != eH.Info[key] {
			t.Errorf("problem with header Info")
		}
	}
	for key := range aH.Filter {
		if aH.Filter[key] != eH.Filter[key] {
			t.Errorf("problem with header Filter")
		}
	}
	for key := range aH.Format {
		if aH.Format[key] != eH.Format[key] {
			t.Errorf("problem with header Format")
		}
	}
	for key := range aH.Samples {
		if aH.Samples[key] != eH.Samples[key] {
			t.Errorf("problem with header chroms")
		}
	}
}
