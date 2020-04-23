package alleles

import (
	"fmt"
	"github.com/vertgenlab/gonomics/fasta"
	"github.com/vertgenlab/gonomics/simpleGraph"
	"sync"
)

func SyncAlleleStreams(reference interface{}, alleleStreams ...<-chan *Allele) <-chan []*Allele {
	fmt.Println("\nSTART SYNC ALLELES")
	answer := make(chan []*Allele)

	var currLoc *Location = &Location{"dummy", 0}

	var receiveAlleles = make(chan *Allele)

	var wg sync.WaitGroup
	locks := make([]*bool, 0)
	locker := sync.NewCond(&sync.Mutex{})

	for _, singleStream := range alleleStreams {
		var lock bool = false
		go sendAlleles(singleStream, receiveAlleles, &lock, locker, &wg, currLoc)
		locks = append(locks, &lock)
	}

	go syncAlleles(reference, receiveAlleles, answer, locks, locker, &wg, currLoc)

	return answer
}

func sendAlleles(inStream <-chan *Allele, outStream chan<- *Allele, lock *bool, locker *sync.Cond, wg *sync.WaitGroup, currLoc *Location) {
	for allele := range inStream {
		// Loop forever until allele is sent on outStream
		for {
			if *lock == true  {
				locker.L.Lock()
				fmt.Println("waiting")
				locker.Wait()
			}
			if allele.Location.Pos == currLoc.Pos && allele.Location.Chr == currLoc.Chr {
				outStream <- allele
				fmt.Println("sent", allele.Location)
				break
			} else {
				// If allele is not ready to send then activate its lock and sleep
				fmt.Println("locked on", allele.Location, currLoc)
				*lock = true
				//wg.Wait()
			}
		}
	}
}

func syncAlleles(reference interface{}, outStream <-chan *Allele, answer chan<- []*Allele, locks []*bool, locker *sync.Cond, wg *sync.WaitGroup, currLoc *Location) {

	var curr []*Allele
	var k int
	var refIdx int = 0

	go func(){
		for allele := range outStream {
			curr = append(curr, allele)
		}
	}()

	//locker.L.Lock()
	for nextPos(reference, currLoc, &refIdx) {
		fmt.Println("current pos", currLoc)
		locker.Broadcast()
		//locker.L.Unlock()
		wg.Add(1)
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
				fmt.Println("resetting locks")
				for k = 0; k < len(locks); k++ {
					*locks[k] = false
				}

				if curr != nil {
					answer <- curr
				}
				curr = nil
				wg.Done()
				break
			} else {locker.Broadcast()}
		}
		//locker.L.Lock()
	}
}

func nextPos (reference interface{}, currLoc *Location, refIdx *int) bool {
	switch r := reference.(type) {
	case []*fasta.Fasta:

		if currLoc.Chr == "dummy" {
			currLoc.Chr = r[0].Name
			currLoc.Pos = 0
			return true
		} else {
			if int(currLoc.Pos) < len(r[*refIdx].Seq) { // Increment Position
				currLoc.Pos++
				return true
			} else if *refIdx < len(r){ // Increment Chromosome
				*refIdx++
				currLoc.Chr = r[*refIdx].Name
				currLoc.Pos = 0
				return true
			} else { // Report end of reference
				return false
			}
		}

	case *simpleGraph.SimpleGraph:

		if currLoc.Chr == "dummy" {
			currLoc.Chr = r.Nodes[0].Name
			currLoc.Pos = 0
			return true
		} else {
			if int(currLoc.Pos) < len(r.Nodes[*refIdx].Seq) { // Increment Position
				currLoc.Pos++
				return true
			} else if *refIdx < len(r.Nodes) { // Increment Node
				*refIdx++
				currLoc.Chr = r.Nodes[*refIdx].Name
				currLoc.Pos = 0
				return true
			} else { // Report end of reference
				return false
			}
		}
	}
	return false
}