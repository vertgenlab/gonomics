package main

import (
	"github.com/vertgenlab/gonomics/exception"
	"github.com/vertgenlab/gonomics/fileio"
	"io/ioutil"
	"os"
	"testing"

	"github.com/vertgenlab/gonomics/bed"
)

var IntervalOverlapTests = []struct {
	InFile       string
	OutFile      string
	SelectFile   string
	NonOverlap   bool
	Threads      int
	Aggregate    bool
	Relationship string
	MergedOutput bool
	ExpectedFile string
}{
	{InFile: "testdata/test.bed",
		OutFile:      "testdata/out.bed",
		SelectFile:   "testdata/test.vcf",
		NonOverlap:   false,
		Threads:      1,
		Aggregate:    false,
		Relationship: "any",
		MergedOutput: false,
		ExpectedFile: "testdata/expected.bed",
	},
	{InFile: "testdata/test.bed",
		OutFile:      "testdata/out.mergedOutput.bed",
		SelectFile:   "testdata/test.vcf",
		NonOverlap:   false,
		Threads:      1,
		Aggregate:    false,
		Relationship: "any",
		MergedOutput: true,
		ExpectedFile: "testdata/expected.mergedOutput.bed",
	},
	{InFile: "testdata/test.bed",
		OutFile:      "testdata/out.nonOverlap.bed",
		SelectFile:   "testdata/test.vcf",
		NonOverlap:   true,
		Threads:      1,
		Aggregate:    false,
		Relationship: "any",
		MergedOutput: false,
		ExpectedFile: "testdata/expected.nonOverlap.bed",
	},
}

func TestIntervalOverlap(t *testing.T) {
	var s *Settings
	var err error
	var out *fileio.EasyWriter
	var answer chan *queryAnswer
	for _, v := range IntervalOverlapTests {
		s = &Settings{
			Input:        v.InFile,
			Output:       v.OutFile,
			SelectFile:   v.SelectFile,
			NonOverlap:   v.NonOverlap,
			Threads:      v.Threads,
			Aggregate:    v.Aggregate,
			Relationship: v.Relationship,
			MergedOutput: v.MergedOutput,
		}
		answer = intervalOverlap(s)
		out = fileio.EasyCreate(v.OutFile)
		writeToFile(answer, out, v.MergedOutput, v.NonOverlap)
		err = out.Close()
		exception.PanicOnErr(err)

		if !fileio.AreEqual(v.OutFile, v.ExpectedFile) {
			t.Errorf("Error: intervalOverlap test output is not as expeected.")
		} else {
			err = os.Remove(v.OutFile)
			exception.PanicOnErr(err)
		}
	}
}

func BenchmarkAssertion(b *testing.B) {
	options := &Settings{
		Input:      "testdata/test.vcf",
		Output:     "/dev/stdout",
		SelectFile: "testdata/test.bed",
		//Extend:          0,
		NonOverlap: false,
		Threads:    1,
		//PercentOverlap:  0,
		//BaseOverlap:     0,
		Aggregate:    false,
		Relationship: "any",
		MergedOutput: false,
		//SwapTargetQuery: false,
	}
	for i := 0; i < b.N; i++ {
		answer := intervalOverlap(options)
		writeToFile(answer, ioutil.Discard, options.MergedOutput, options.NonOverlap)
	}
}

func BenchmarkAssertWrite(b *testing.B) {
	for i := 0; i < b.N; i++ {
		b.StopTimer()
		a := bed.Bed{FieldsInitialized: 3}
		var i interface{}
		i = a
		b.StartTimer()
		i.(fileWriter).WriteToFileHandle(ioutil.Discard)
	}
}

func BenchmarkAllocWrite(b *testing.B) {
	for i := 0; i < b.N; i++ {
		b.StopTimer()
		a := bed.Bed{FieldsInitialized: 3}
		var i fileWriter
		i = a
		b.StartTimer()
		i.WriteToFileHandle(ioutil.Discard)
	}
}
