package sam

import (
	"github.com/vertgenlab/gonomics/fileio"
	"log"
	"strings"
	"testing"
)

func Test10xBarcode(t *testing.T) {
	TenXPrettyPrint("testdata/tenXbarcodeTest.sam")
}

func TestReadingBarcode(t *testing.T) {
	log.SetFlags(log.Ltime)
	reader := make(chan *SamAln)
	samfile := fileio.EasyOpen("testdata/tenXbarcodeTest.sam")
	defer samfile.Close()
	ReadHeader(samfile)
	go ReadToChan(samfile, reader)
	var bxTag string
	for read := range reader {
		bxTag = LinkedReadsBarcode(read)
		if strings.HasPrefix(bxTag, "BX") || strings.HasPrefix(bxTag, "RX") {
			log.Printf("Pass: %s\n", bxTag)
		} else {
			t.Errorf("Error: Get 10x barcode functions not working properly or not barcode exists\n")
		}
	}
}
