package sam

import (
	"github.com/vertgenlab/gonomics/giraf"
	"strings"
	//"github.com/vertgenlab/gonomics/dna"
	"fmt"
	"github.com/vertgenlab/gonomics/cigar"
	"github.com/vertgenlab/gonomics/fileio"
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

func LinkedReadsBarcode(samfile *SamAln) string {
	barcodeInfo := BarcodeInfoMap(samfile)
	_, ok := barcodeInfo["BX"]
	if ok {
		return customNoteToString(barcodeInfo["BX"])
	} else {
		return customNoteToString(barcodeInfo["RX"])
	}
}

func customNoteToString(n *giraf.Note) string {
	return fmt.Sprintf("%s_%s", n.Tag, n.Value[:len(n.Value)-2])
}

func switchBxTag(samfile *SamAln) *SamAln {
	//bxTag := LinkedReadsBarcode(samfile)
	//read := samfile.QName
	samfile.Extra = fmt.Sprintf("RX:Z:%s\t%s", samfile.QName, samfile.Extra)
	samfile.QName = LinkedReadsBarcode(samfile) //fmt.Sprintf("%s", bxTag)
	//infoMap := BarcodeInfoMap(samfile)
	//tags := fmt.Sprintf("\tXR:Z:%s", samfile.RName)
	//for keys := range infoMap {
	//	if keys == "BX" {
	//		bxTag := fmt.Sprintf("%s_%s", infoMap[keys].Tag, infoMap[keys].Value[:len(infoMap[keys].Value)-2])
	//samfile.RName, samfile.RNext  = bxTag, bxTag
	//	}
	//}
	//samfile.Extra = tags
	return samfile

}

func bxTagWorker(reader <-chan *SamAln, writer chan<- *SamAln, wg *sync.WaitGroup) {
	for align := range reader {
		writer <- switchBxTag(align)
	}
	wg.Done()
}

func SwitchBxTagChannel(filename string, output string, threads int) {
	samfile := fileio.EasyOpen(filename)
	defer samfile.Close()
	header := ReadHeader(samfile)

	var worker, toFile sync.WaitGroup
	reader, writer := make(chan *SamAln), make(chan *SamAln)
	go ReadToChan(samfile, reader)
	worker.Add(threads)
	for i := 0; i < threads; i++ {
		go bxTagWorker(reader, writer, &worker)
	}
	toFile.Add(1)
	go SamChanToFile(writer, output, header, &toFile)

	worker.Wait()
	close(writer)
	toFile.Wait()

	//	modified := fileio.EasyCreate(output)
	//	WriteHeaderToFileHandle(modified, header)
	//var tenX *SamAln
	//	for read, done := NextAlignment(samFile); done != true; read, done = NextAlignment(samFile) {
	//		tenX := switchBxTag(read)
	//		WriteAlnToFileHandle(modified, tenX)
	//	}

}

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

func getHeader(filename string) *SamHeader {
	file := fileio.EasyOpen(filename)
	header := ReadHeader(file)
	defer file.Close()
	return header
}

func TenXPrettyPrint(filename string) {
	samfile := fileio.EasyOpen(filename)
	defer samfile.Close()
	ReadHeader(samfile)
	linkedReads := make(chan *SamAln)
	go ReadToChan(samfile, linkedReads)
	for barcodes := range linkedReads {
		if !IsForwardRead(barcodes) {
			fmt.Printf("%s\t%s\t%s\t%s\n", barcodes.RName, barcodes.QName, LinkedReadsBarcode(barcodes), cigar.ToString(barcodes.Cigar))
			//fmt.Printf("%s\n", )
		} else {

		}
	}
}
