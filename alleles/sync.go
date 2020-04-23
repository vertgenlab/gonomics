package alleles

import (
	"github.com/vertgenlab/gonomics/fasta"
	"github.com/vertgenlab/gonomics/simpleGraph"
	"sync"
)

func SyncAlleleStreams(reference interface{}, alleleStreams ...<-chan *Allele) <-chan []*Allele {
	answer := make(chan []*Allele)

	var currLoc *Location //TODO: how to define start location???
	var receiveAlleles = make(chan *Allele)

	var wg sync.WaitGroup
	locks := make([]*bool, 0)

	for _, singleStream := range alleleStreams {
		var lock bool = false
		go sendAlleles(singleStream, receiveAlleles, &lock, &wg, currLoc)
		locks = append(locks, &lock)
	}

	go syncAlleles(reference, receiveAlleles, answer, locks, &wg, currLoc)

	return answer
}

func sendAlleles(inStream <-chan *Allele, outStream chan<- *Allele, lock *bool, wg *sync.WaitGroup, currLoc *Location) {
	for allele := range inStream {
		// Loop forever until allele is sent on outStream
		for {
			if allele.Location.Pos == currLoc.Pos && allele.Location.Chr == currLoc.Chr {
				outStream <- allele
				break
			} else {
				// If allele is not ready to send then activate its lock and sleep
				*lock = true
				wg.Wait()
			}
		}
	}
}

func nextPos (reference interface{}, currLoc *Location) bool {
	switch reference.(type) {
	case []*fasta.Fasta:

		// if end of file return false
	case *simpleGraph.SimpleGraph:

		// if end of file return false
	}
	return false
}

func syncAlleles(reference interface{}, outStream <-chan *Allele, answer chan<- []*Allele, locks []*bool, wg *sync.WaitGroup, currLoc *Location) {

	var curr []*Allele
	var k int

	go func(){
		for allele := range outStream {
			curr = append(curr, allele)
		}
	}()

	for nextPos(reference, currLoc) {
		for {
			var readyToProgress bool = true
			// Check to see if all workers are ready for next currLocation
			for k = 0; k < len(locks); k++ {
				if *locks[k] == false {
					readyToProgress = false
					break
				}
			}
			// If all workers are ready, then reset the locks, send the allele, and progress
			if readyToProgress == true {
				for k = 0; k < len(locks); k++ {
					*locks[k] = false
				}

				if curr != nil {
					answer <- curr
				}
				curr = nil
				wg.Done()
				break
			}
		}
	}
}