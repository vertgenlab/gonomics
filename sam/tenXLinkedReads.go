package sam

import (
	"fmt"
	"github.com/vertgenlab/gonomics/cigar"
	"github.com/vertgenlab/gonomics/giraf"
	"strings"
	"sync"
)

func BarcodeInfoMap(samfile *SamAln) map[string]*giraf.Note {
	words := strings.Split(samfile.Extra, "\t")
	answer := make(map[string]*giraf.Note)
	var values []string
	for _, tag := range words {
		values = strings.SplitN(tag, ":", 3)
		answer[values[0]] = &giraf.Note{Tag: values[0], Type: []rune(values[1])[0], Value: values[2]}
	}
	return answer
}

//Barcodes corrected for sequence errors are taged BX, while uncorrected barcodes are labeled RX
//BX tags do not exist for all reads
func LinkedReadsBarcode(samLine *SamAln) string {
	barcodeInfo := BarcodeInfoMap(samLine)
	data, ok := barcodeInfo["BX"]
	if ok {
		return customNoteToString(data)
	} else {
		return customNoteToString(barcodeInfo["RX"])
	}
}

func customNoteToString(n *giraf.Note) string {
	return fmt.Sprintf("%s_%s", n.Tag, n.Value[:len(n.Value)-2])
}

func switchBxTag(samLine *SamAln) *SamAln {
	samLine.Extra = fmt.Sprintf("RX:Z:%s\t%s", samLine.QName, samLine.Extra)
	samLine.QName = LinkedReadsBarcode(samLine) //fmt.Sprintf("%s", bxTag)
	return samLine

}

func BxTagWorker(reader <-chan *SamAln, writer chan<- *SamAln, wg *sync.WaitGroup) {
	for align := range reader {
		writer <- switchBxTag(align)
	}
	wg.Done()
}

//TODO modify a header containing BX barcodes as "reference"
/*
func BarcodeHeader(filename string, wg *sync.WaitGroup) *SamHeader {
	//wg.Add(1)
	reader := make(chan *SamAln)
	var header SamHeader
	header.Text = append(header.Text, "@HD\tVN:1.6\tSO:unsorted")
	go ReadToChan(filename, reader)
	bxMap := make(map[string]string)
	var key string
	for align := range reader {
		key = LinkedReadsBarcode(align)
		_, ok := bxMap[key]
		if !ok {
			bxMap[key] = key
		}
	}
	for tenX := range bxMap {
		header.Text = append(header.Text, fmt.Sprintf("@SQ\tSN:%s\tLN:%d", tenX, 16))
	}
	wg.Done()
	return &header
}*/

func TenXPrettyPrint(filename string) {
	linkedReads, _ := GoReadToChan(filename)
	for barcodes := range linkedReads {
		if !IsForwardRead(barcodes) {
			fmt.Printf("%s\t%s\t%s\t%s\n", barcodes.RName, barcodes.QName, LinkedReadsBarcode(barcodes), cigar.ToString(barcodes.Cigar))
		}
	}
}
