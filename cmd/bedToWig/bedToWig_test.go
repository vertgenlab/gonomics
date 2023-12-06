package main

import (
	"os"
	"testing"

	"github.com/vertgenlab/gonomics/exception"
	"github.com/vertgenlab/gonomics/fileio"
)

var bedToWigTests = []struct {
	InFile          string
	RefFile         string
	OutFile         string
	ExpectedFile    string
	Method          string
	Missing         float64
	UseRange        bool
	AnnotationField int
}{
	{InFile: "testdata/test.bed",
		RefFile:         "testdata/ref.chrom.sizes",
		OutFile:         "testdata/test.Score.wig",
		ExpectedFile:    "testdata/score.Expected.wig",
		Method:          "Score",
		Missing:         0,
		UseRange:        false,
		AnnotationField: 0},
	{InFile: "testdata/test.bed",
		RefFile:         "testdata/ref.chrom.sizes",
		OutFile:         "testdata/test.Reads.wig",
		ExpectedFile:    "testdata/reads.Expected.wig",
		Method:          "Reads",
		Missing:         0,
		UseRange:        false,
		AnnotationField: 0},
	{InFile: "testdata/test.bed",
		RefFile:         "testdata/ref.chrom.sizes",
		OutFile:         "testdata/test.Name.wig",
		ExpectedFile:    "testdata/name.Expected.wig",
		Method:          "Name",
		Missing:         0,
		UseRange:        false,
		AnnotationField: 0},
	{InFile: "testdata/test.bed",
		RefFile:         "testdata/ref.chrom.sizes",
		OutFile:         "testdata/test.missing.Name.wig",
		ExpectedFile:    "testdata/name.missing.Expected.wig",
		Method:          "Name",
		Missing:         -1.0,
		UseRange:        false,
		AnnotationField: 0},
	{InFile: "testdata/test.range.bed",
		RefFile:         "testdata/ref.chrom.sizes",
		OutFile:         "testdata/test.range.Name.wig",
		ExpectedFile:    "testdata/name.range.Expected.wig",
		Method:          "Name",
		Missing:         -1.0,
		UseRange:        true,
		AnnotationField: 0},
	{InFile: "testdata/test.range.bed",
		RefFile:         "testdata/ref.chrom.sizes",
		OutFile:         "testdata/test.range.Score.wig",
		ExpectedFile:    "testdata/score.range.Expected.wig",
		Method:          "Score",
		Missing:         -1.0,
		UseRange:        true,
		AnnotationField: 0},
	{InFile: "testdata/test.annotation.bed",
		RefFile:         "testdata/annotation.chrom.sizes",
		OutFile:         "testdata/test.Annotation.wig",
		ExpectedFile:    "testdata/expected.Annotation.wig",
		Method:          "Annotation",
		Missing:         -1.0,
		UseRange:        false,
		AnnotationField: 0},
	{InFile: "testdata/test.annotation.bed",
		RefFile:         "testdata/annotation.chrom.sizes",
		OutFile:         "testdata/test.Annotation.Field2.wig",
		ExpectedFile:    "testdata/expected.Annotation.Field2.wig",
		Method:          "Annotation",
		Missing:         -1.0,
		UseRange:        false,
		AnnotationField: 2},
}

func TestBedToWig(t *testing.T) {
	var err error
	var s Settings
	for _, v := range bedToWigTests {
		s = Settings{
			Method:          v.Method,
			InFile:          v.InFile,
			RefFile:         v.RefFile,
			OutFile:         v.OutFile,
			Missing:         v.Missing,
			UseRange:        v.UseRange,
			AnnotationField: v.AnnotationField,
		}
		bedToWig(s)
		if !fileio.AreEqual(v.OutFile, v.ExpectedFile) {
			t.Errorf("Error: output of bedToWig was not as expected.")
		} else {
			err = os.Remove(v.OutFile)
			exception.PanicOnErr(err)
		}
	}
}
