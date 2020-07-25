package bam

import (
	"github.com/vertgenlab/gonomics/sam"
	"testing"
	//"log"
)

func TestBgzipBlocks(t *testing.T) {

	bamFile := Read("testdata/tenXbarcodeTest.bam")
	samFile, _ := sam.Read("testdata/tenXbarcodeTest.sam")
	if len(bamFile) != len(samFile.Aln) {
		t.Errorf("Error: File lines are not equal...\n")
	}
	for i := 0; i < len(bamFile); i++ {
		if !sam.IsEqualDebug(bamFile[i], samFile.Aln[i]) {
			t.Fatalf("Error: Did not create the same sam file as samtools view...\n")
		}
	}
	//Read("testdata/small.bam")
	//bam := NewBamReader("testdata/rabsCanuHicTenX.gasAcu1.bam")
	//defer bam.File.Close()
	//ReadHeader(bam)
	//bamTryAgain(bam)
	//bamLine(bam)
	//ans := processBamRecord(bam)
	//bamBlocks(ans)
	//reader := GunzipReader(file.File)
	//ReadHeader(file.File)
	//log.Printf("%s\n", reader.Header.Text)
}
