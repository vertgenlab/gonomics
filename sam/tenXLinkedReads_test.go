package sam

import (
	"log"
	"strings"
	"testing"
) /*
//TODO: Only debugging prints, need to come up with a better test function
func Test10xBarcode(t *testing.T) {
	TenXPrettyPrint("testdata/tenXbarcodeTest.sam")
}*/
//TODO hard code this into a small toy example
func TestReadingBarcode(t *testing.T) {
	log.SetFlags(log.Ltime)
	reader, _ := GoReadToChan("testdata/tenXbarcodeTest.sam")
	var bxTag string
	for read := range reader {
		bxTag = LinkedReadsBarcode(read)
		if strings.HasPrefix(bxTag, "BX") || strings.HasPrefix(bxTag, "RX") {
			//log.Printf("Pass: %s\n", bxTag)
		} else {
			t.Errorf("Error: Get 10x barcode functions not working properly or not barcode exists\n")
		}
	}
}
