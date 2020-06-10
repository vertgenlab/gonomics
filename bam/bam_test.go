package bam

import (
	"testing"
	//"github.com/vertgenlab/gonomics/fileio"
	//"log"

)

func TestBgzipBlocks(t *testing.T) {
	bam := NewBamReader("testdata/rabsToGasAcu1.bam")
	ReadHeader(bam)
	processBamRecord(bam)
	//bamBlocks(ans)
	//reader := GunzipReader(file.File)
	//ReadHeader(file.File)
	//log.Printf("%s\n", reader.Header.Text)
}