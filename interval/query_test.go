package interval

import (
	"bytes"
	"log"
	"os"
	"strings"
	"testing"

	"github.com/vertgenlab/gonomics/bed"
	"github.com/vertgenlab/gonomics/fileio"
)

func TestReadToChanBed(t *testing.T) {
	var b1 bed.Bed = bed.Bed{Chrom: "chr1", ChromStart: 100, ChromEnd: 200, Name: "First", Score: 1, Strand: '+', FieldsInitialized: 6}
	var beds []bed.Bed = []bed.Bed{b1}
	var bedFile string = "testdata/interval.bed"

	bed.Write(bedFile, beds)

	// Create the interval channel
	interval := make(chan Interval, 1)
	ReadToChan(bedFile, interval)

	// Check the intervals read from the channel
	expectedBeds := []bed.Bed{b1}
	for _, expected := range expectedBeds {
		val, ok := <-interval
		if !ok {
			t.Fatalf("Expected more intervals, but channel was closed")
		}
		bedVal, ok := val.(*bed.Bed)
		if !ok {
			t.Fatalf("Expected *bed.Bed type, got %T", val)
		}
		if bedVal.Chrom != expected.Chrom || bedVal.ChromStart != expected.ChromStart || bedVal.ChromEnd != expected.ChromEnd {
			t.Errorf("Error: %v != %v", bedVal, expected)
		}
	}

	// Verify that no more intervals are left in the channel
	if _, ok := <-interval; ok {
		t.Errorf("Expected channel to be closed, but it still has values")
	}
	fileio.EasyRemove(bedFile)
}

func TestReadToChanUnknownFileType(t *testing.T) {
	if os.Getenv("TEST_FATAL") != "1" {
		return
	}
	// Create a temporary file with an unknown extension
	unknownData := "testdata/interval.unknown"
	unknownFile := fileio.EasyCreate(unknownData)
	defer unknownFile.Close()
	defer fileio.EasyRemove(unknownData)

	// Channel to receive intervals (not expected to be used)
	intervalCh := make(chan Interval)

	// Call ReadToChan in a goroutine
	ReadToChan(unknownData, intervalCh)
	// Capture log output
	var buf bytes.Buffer
	log.SetOutput(&buf)
	defer log.SetOutput(os.Stderr)
	expectedMessage := "Error:"
	logOutput := buf.String()
	if !strings.Contains(expectedMessage, logOutput) {
		t.Errorf("Expected log message:\n%s\nNot found in actual log:\n%s", expectedMessage, logOutput)
	} else {
		t.Log("log.Fatalf called with expected message - Test Passed!")
	}
}
