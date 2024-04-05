package main

import (
	"fmt"
	"github.com/vertgenlab/gonomics/exception"
	"github.com/vertgenlab/gonomics/fileio"
	"github.com/vertgenlab/gonomics/numbers/parse"
	"strings"
	"sync"
)

type multiReads struct {
	qname string
	key   string
	pVal  float64
}

type uniReadBin struct {
	key   string
	count int
}

func goReadUniReadBin(inFile string) <-chan uniReadBin {
	var wg sync.WaitGroup
	data := make(chan uniReadBin, 1000)
	wg.Add(1)
	go readUniReadBinToChan(inFile, data, &wg)

	go func() {
		wg.Wait()
		close(data)
	}()

	return data
}

func readUniReadBinToChan(inFile string, data chan<- uniReadBin, wg *sync.WaitGroup) {
	var line uniReadBin
	var done bool
	in := fileio.EasyOpen(inFile)

	for line, done = nextUniReadBin(in); !done; line, done = nextUniReadBin(in) {
		data <- line
	}
	err := in.Close()
	exception.PanicOnErr(err)
	wg.Done()
}

func nextUniReadBin(reader *fileio.EasyReader) (uniReadBin, bool) {
	line, done := fileio.EasyNextRealLine(reader)
	if done {
		return uniReadBin{}, true
	}
	return processsUniReadBin(line), false
}

func processsUniReadBin(line string) uniReadBin {
	var u uniReadBin
	slc := strings.Split(line, "\t")
	u.count = parse.StringToInt(slc[4])
	u.key = strings.Join(slc[0:4], "_")
	return u
}

func goReadMultiReadToChan(inFile string) <-chan multiReads {
	var wg sync.WaitGroup
	data := make(chan multiReads, 1000)
	wg.Add(1)
	go readMultiReadToChan(inFile, data, &wg)

	go func() {
		wg.Wait()
		close(data)
	}()

	return data
}

func readMultiReadToChan(inFile string, data chan<- multiReads, wg *sync.WaitGroup) {
	var line multiReads
	var done bool
	in := fileio.EasyOpen(inFile)

	for line, done = nextMultiRead(in); !done; line, done = nextMultiRead(in) {
		data <- line
	}
	err := in.Close()
	exception.PanicOnErr(err)
	wg.Done()
}

func nextMultiRead(reader *fileio.EasyReader) (multiReads, bool) {
	line, done := fileio.EasyNextRealLine(reader)
	if done {
		return multiReads{}, true
	}
	return processsMultiReads(line), false
}

func processsMultiReads(line string) multiReads {
	var mr multiReads
	slc := strings.Split(line, "\t")
	mr.qname = slc[0]
	mr.key = strings.Join(slc[1:5], "_")
	mr.pVal = parse.StringToFloat64(slc[5])
	return mr
}

func addAssignedReads(uniReadBins, assignedMultiReads, outFile string) {
	mp := makeMultiReadMap(assignedMultiReads)
	updateUniReadBins(uniReadBins, mp, outFile)
}

func updateUniReadBins(uniReadBins string, mp map[string]int, outFile string) {
	o := fileio.EasyCreate(outFile)
	ch := goReadUniReadBin(uniReadBins)
	for i := range ch {
		i.count += mp[i.key]
		writeUpdatedCount(i, o)
	}
	err := o.Close()
	exception.PanicOnErr(err)
}

func writeUpdatedCount(u uniReadBin, o *fileio.EasyWriter) {
	slc := strings.Split(u.key, "_")
	fileio.WriteToFileHandle(o, fmt.Sprintf("%s\t%s\t%s\t%s\t%d", slc[0], slc[1], slc[2], slc[3], u.count))
}

func makeMultiReadMap(multiReads string) map[string]int {
	mp := make(map[string]int)
	mr := goReadMultiReadToChan(multiReads)

	for i := range mr {
		if i.pVal <= 0.5 {
			continue
		}
		mp[i.key]++
	}
	return mp
}

func main() {
	uniReads := "/Users/sethweaver/Downloads/gonomicsMHIC/testdata/output/uniContacts.txt"
	multiReadList := "/Users/sethweaver/Downloads/gonomicsMHIC/testdata/output/testOutS6.go.txt"
	outFile := "/Users/sethweaver/Downloads/gonomicsMHIC/testdata/output/uniMultiContactBins.txt"

	addAssignedReads(uniReads, multiReadList, outFile)
}
