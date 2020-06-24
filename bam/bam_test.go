package bam

import (
	"testing"
	//"github.com/vertgenlab/gonomics/fileio"
	//"log"

)

func TestBgzipBlocks(t *testing.T) {
	bam := NewBamReader("testdata/rabsCanuHicTenX.gasAcu1.bam")
	defer bam.File.Close()
	ReadHeader(bam)
	bamLine(bam)
	//ans := processBamRecord(bam)
	//bamBlocks(ans)
	//reader := GunzipReader(file.File)
	//ReadHeader(file.File)
	//log.Printf("%s\n", reader.Header.Text)
}