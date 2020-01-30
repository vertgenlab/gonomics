package routines

import (
	"github.com/vertgenlab/gonomics/fileio"
	"github.com/vertgenlab/gonomics/sam"
	"log"
	"sync"
)

type SamOperation func(*sam.SamAln, interface{}) interface{}

func goReadSam(samFilename string, sendLine chan *sam.SamAln, fileType string) {
	samFile := fileio.EasyOpen(samFilename)
	defer samFile.Close()
	var done = false
	var aln *sam.SamAln



	sam.ReadHeader(samFile)

	progressMeter := 0
	for aln, done = sam.NextAlignment(samFile); done != true; aln, done = sam.NextAlignment(samFile) {

		progressMeter++
		if progressMeter % 10000 == 0 {
			log.Println("Processed", progressMeter, "Lines")
		}
		sendLine <- aln
	}
	close(sendLine)
}

func ReadByLine(samFilename string, function SamOperation, data1 interface{}, threads int, fileType string) chan interface{}{
	var wg sync.WaitGroup
	records := make(chan *sam.SamAln)
	output := make(chan interface{})

	go goReadSam(samFilename, records, fileType)

	for k := 0; k < threads; k++ {
		wg.Add(1)
		go func() {
			for data := range records {
				a := function(data, data1)
				output <- a
			}
			wg.Done()
			return
		}()
	}

	go func() {
		wg.Wait()
		close(output)
	}()

	return output
}



func ReadPos(read *sam.SamAln, data1 interface{}) interface{} {
	var data = data1.(*input)
	data.a = 9
	return read.Seq
}