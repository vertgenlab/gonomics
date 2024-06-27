package interval

import (
	"bytes"
	"fmt"
	"log"
	"os"
	"strings"
	"testing"

	"github.com/vertgenlab/gonomics/axt"
	"github.com/vertgenlab/gonomics/bed"
	"github.com/vertgenlab/gonomics/chain"
	"github.com/vertgenlab/gonomics/cigar"
	"github.com/vertgenlab/gonomics/dna"
	"github.com/vertgenlab/gonomics/fileio"
	"github.com/vertgenlab/gonomics/sam"
	"github.com/vertgenlab/gonomics/vcf"
)

func TestReadToChanAxt(t *testing.T) {
	var a axt.Axt = axt.Axt{
		RName:      "TargetGenome",
		RStart:     1,
		REnd:       10,
		QName:      "QueryGenome",
		QStart:     1,
		QEnd:       10,
		QStrandPos: false,
		RSeq:       dna.StringToBases("CGTCCCTCGA"),
		QSeq:       dna.StringToBases("CGTCCCTCGA"),
	}
	var axts []axt.Axt = []axt.Axt{a}
	var axtFile string = "testdata/interval.axt"

	axt.Write(axtFile, axts)

	// Create the interval channel
	interval := GoReadToChan(axtFile)

	// Check the intervals read from the channel
	expected := []axt.Axt{a}
	for _, expected := range expected {
		val, ok := <-interval
		if !ok {
			t.Fatalf("Expected more intervals, but channel was closed")
		}
		axtVal, ok := val.(axt.Axt)
		if !ok {
			t.Fatalf("Expected *axt.Axt type, got %T", val)
		}
		if axtVal.GetChrom() != expected.GetChrom() || axtVal.GetChromStart() != expected.GetChromStart() || axtVal.GetChromEnd() != expected.GetChromEnd() {
			t.Errorf("Error: %v != %v", axtVal, expected)
		}
	}

	// Verify that no more intervals are left in the channel
	if _, ok := <-interval; ok {
		t.Errorf("Expected channel to be closed, but it still has values")
	}
	fileio.EasyRemove(axtFile)
}

func TestReadToChanBed(t *testing.T) {
	var b bed.Bed = bed.Bed{Chrom: "chr1", ChromStart: 100, ChromEnd: 200, Name: "First", Score: 1, Strand: '+', FieldsInitialized: 6}
	var beds []bed.Bed = []bed.Bed{b}
	var bedFile string = "testdata/interval.bed"

	bed.Write(bedFile, beds)

	// Create the interval channel
	interval := GoReadToChan(bedFile)

	// Check the intervals read from the channel
	expectedBeds := []bed.Bed{b}
	for _, expected := range expectedBeds {
		val, ok := <-interval
		if !ok {
			t.Fatalf("Expected more intervals, but channel was closed")
		}
		bedVal, ok := val.(*bed.Bed)
		if !ok {
			t.Fatalf("Expected *bed.Bed type, got %T", val)
		}
		if bedVal.GetChrom() != expected.GetChrom() || bedVal.GetChromStart() != expected.GetChromStart() || bedVal.GetChromEnd() != expected.GetChromEnd() {
			t.Errorf("Error: %v != %v", bedVal, expected)
		}
	}

	// Verify that no more intervals are left in the channel
	if _, ok := <-interval; ok {
		t.Errorf("Expected channel to be closed, but it still has values")
	}
	fileio.EasyRemove(bedFile)
}

func TestReadToChanChain(t *testing.T) {
	var c chain.Chain = chain.Chain{
		Score:     4766,
		TName:     "chrI",
		TSize:     600,
		TStrand:   true,
		TStart:    550,
		TEnd:      600,
		QName:     "contig_12",
		QSize:     50,
		QStrand:   false,
		QStart:    0,
		QEnd:      50,
		Alignment: make([]chain.BaseStats, 1),
		Id:        1,
	}

	var chains []chain.Chain = []chain.Chain{c}
	var chainFile string = "testdata/interval.chain"
	chain.Write(chainFile, chains, chain.HeaderComments{})

	// Create the interval channel
	interval := make(chan Interval, 1)
	ReadToChan(chainFile, interval)

	// Check the intervals read from the channel
	expected := []chain.Chain{c}
	for _, expected := range expected {
		val, ok := <-interval
		if !ok {
			t.Fatalf("Expected more intervals, but channel was closed")
		}
		axtVal, ok := val.(chain.Chain)
		if !ok {
			t.Fatalf("Expected chain.Chain type, got %T", val)
		}
		if axtVal.GetChrom() != expected.GetChrom() || axtVal.GetChromStart() != expected.GetChromStart() || axtVal.GetChromEnd() != expected.GetChromEnd() {
			t.Errorf("Error: %v != %v", axtVal, expected)
		}
	}
	// Verify that no more intervals are left in the channel
	if _, ok := <-interval; ok {
		t.Errorf("Expected channel to be closed, but it still has values")
	}
	fileio.EasyRemove(chainFile)
}

func TestReadToChanSam(t *testing.T) {
	var s sam.Sam = sam.Sam{
		QName: "r001",
		Flag:  99,
		MapQ:  30,
		RName: "ref",
		Pos:   7,
		Cigar: cigar.FromString("8M2I4M1D3M"),
		RNext: "=",
		PNext: 37,
		TLen:  39,
		Seq:   dna.StringToBases("TTAGATAAAGGATACTG"),
		Qual:  "*",
		Extra: "",
	}

	var sams []sam.Sam = []sam.Sam{s}
	var samFile string = "testdata/interval.sam"
	sam.Write(samFile, sams, sam.Header{})

	// Create the interval channel
	interval := make(chan Interval, 1)
	ReadToChan(samFile, interval)

	// Check the intervals read from the channel
	expected := []sam.Sam{s}
	for _, expected := range expected {
		val, ok := <-interval
		if !ok {
			t.Fatalf("Expected more intervals, but channel was closed")
		}
		axtVal, ok := val.(sam.Sam)
		if !ok {
			t.Fatalf("Expected sam.Sam type, got %T", val)
		}
		if axtVal.GetChrom() != expected.GetChrom() || axtVal.GetChromStart() != expected.GetChromStart() || axtVal.GetChromEnd() != expected.GetChromEnd() {
			t.Errorf("Error: %v != %v", axtVal, expected)
		}
	}
	// Verify that no more intervals are left in the channel
	if _, ok := <-interval; ok {
		t.Errorf("Expected channel to be closed, but it still has values")
	}
	fileio.EasyRemove(samFile)
}

func TestReadToChanVcf(t *testing.T) {
	var v vcf.Vcf = vcf.Vcf{Chr: "chr1", Pos: 2, Ref: "A", Alt: []string{"T"}}

	var vcfs []vcf.Vcf = []vcf.Vcf{v}
	var vcfFile string = "testdata/interval.vcf.gz"

	vcf.Write(vcfFile, vcfs)

	// Create the interval channel
	interval := make(chan Interval, 1)
	ReadToChan(vcfFile, interval)

	// Check the intervals read from the channel
	expected := []vcf.Vcf{v}
	for _, expected := range expected {
		val, ok := <-interval
		if !ok {
			t.Fatalf("Expected more intervals, but channel was closed")
		}
		axtVal, ok := val.(vcf.Vcf)
		if !ok {
			t.Fatalf("Expected vcf.Vcf type, got %T", val)
		}
		if axtVal.GetChrom() != expected.GetChrom() || axtVal.GetChromStart() != expected.GetChromStart() || axtVal.GetChromEnd() != expected.GetChromEnd() {
			t.Errorf("Error: %v != %v", axtVal, expected)
		}
	}

	// Verify that no more intervals are left in the channel
	if _, ok := <-interval; ok {
		t.Errorf("Expected channel to be closed, but it still has values")
	}
	fileio.EasyRemove(vcfFile)
}

func TestReadToChanUnknownFileType(t *testing.T) {
	defer func() {
		if r := recover(); r == nil {
			t.Errorf("Error: Unknown file type should throw error...\n")
		}
	}()
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

	expectedMessage := fmt.Sprintf("file type of %s is unknown. Does not match any of the following: .bed/.axt/.vcf/sam/.bam/.chain", unknownData)
	logOutput := buf.String()

	if !strings.Contains(expectedMessage, logOutput) {
		t.Errorf("Expected log message:\n%s\nNot found in actual log:\n%s", expectedMessage, logOutput)
	} else {
		t.Log("log.Fatalf called with expected message - Test Passed!")
	}
}
