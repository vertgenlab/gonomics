package main

import (
	"github.com/vertgenlab/gonomics/sam"
	"log"
	"os"
	"sync"
)

func samFilter(filename string, outfile string) {
	sendChan := make(chan sam.Sam, 1000)
	var wg sync.WaitGroup
	reads, header := sam.GoReadToChan(filename)
	go sam.WriteFromChan(sendChan, outfile, header, &wg)
	var totalPairs int
	var pairsRemoved int
	readMap := make(map[string]sam.Sam)
	for s := range reads {
		if _, ok := readMap[s.QName]; !ok {
			totalPairs++
			readMap[s.QName] = s
		} else {
			if sam.MateIsPosStrand(s) && s.Pos >= readMap[s.QName].Pos {
				sendChan <- readMap[s.QName]
				sendChan <- s
			} else if !sam.MateIsPosStrand(s) && s.Pos <= readMap[s.QName].Pos {
				sendChan <- readMap[s.QName]
				sendChan <- s
			} else {
				pairsRemoved++
			}
			delete(readMap, s.QName)
		}
	}
	close(sendChan)
	log.Printf("Removed %d of %d read pairs", pairsRemoved, totalPairs)
}

func main() {
	samFilter(os.Args[1], os.Args[2])
}
