package simulate

import (
	"fmt"
	"github.com/vertgenlab/gonomics/fasta"
	"github.com/vertgenlab/gonomics/dna"
	"github.com/vertgenlab/gonomics/numbers"
	"strconv"
	"math/rand"
)

// GoSimulateAFS takes in a
func GoSimulateAFS (popSize int, mutR float64, genTime int, genomeSize int) []fasta.Fasta {
	allFasta := make([]fasta.Fasta, popSize)
	aPointer := &allFasta
	aSeqPointer := make([]*[]dna.Base, popSize)

	initialSeq := make([]dna.Base, genomeSize)
	var bases = []byte{'A', 'C', 'G', 'T'}
	initialSeq := simulate.RandIntergenicSeq(0.5, genomeSize)
	// Randomly generate the initial genome sequence
	for i := 0; i < genomeSize; i++ {
		initialSeq[i] = dna.ByteToBase(bases[numbers.RandIntInRange(0,4)])
	}

	// Copy the same fasta for the entire population
	for j := 0; j < popSize; j++ {
		(*aPointer)[j].Name = strconv.Itoa(j)
		(*aPointer)[j].Seq = initialSeq
		aSeqPointer[j] = &((*aPointer)[j].Seq)
	}

	// 1st loop through every generation
	for t := 1; t < genTime; t++ {
		var nextFasta = make([]fasta.Fasta, popSize)
		nPointer := &nextFasta
		nSeqPointer := make([]*[]dna.Base, popSize)
		//fmt.Println(allFasta)
		// Empty nextFasta with assigned name
		for x := 0; x < popSize; x++ {
			(*nPointer)[x].Name = strconv.Itoa(x)
			(*nPointer)[x].Seq = make([]dna.Base, genomeSize)
			nSeqPointer[x] = &((*nPointer)[x].Seq)
		}
		
		// 2nd loop through every base 
		for b := 0; b < genomeSize; b++ {
			//fmt.Println(nextFasta)

			// 3rd loop through every individual
			for p := 0; p < popSize; p++ {
				// Which chromosome from previous generation that pth individual inherits from?
				r := numbers.RandIntInRange(0,popSize)
				*(nSeqPointer[p]) = allFasta[r].Seq //randomize index of previous gen
				fmt.Println(b, r, p)
				fmt.Println(allFasta)
				fmt.Println(nextFasta)


				//fmt.Println(r, "---------", p)
				//fmt.Println((*aPointer)[r].Seq )
				//fmt.Println(allFasta[r].Seq)
				
				// If random float generated is less than mutation rate (mutR), mutate the base
				if (rand.Float64() < mutR) {
					//fmt.Println(rand.Float64())
					var cur rune = dna.BaseToRune((*nSeqPointer[p])[b]) // Current Base
					fmt.Println(nextFasta[p].Seq[b])
					fmt.Println((*nSeqPointer[p])[b] == (*aSeqPointer[p])[b])
					switch cur {
					case 'A':
						(*nSeqPointer[p])[b] =  dna.ByteToBase(bases[numbers.RandIntInRange(1,4)])
					case 'C':
						var r int = numbers.RandIntInRange(0,3)
						if r == 1 {
							(*nSeqPointer[p])[b] =  dna.ByteToBase(bases[3])
						} else {
							(*nSeqPointer[p])[b] =  dna.ByteToBase(bases[r])
						}
					case 'G':
						var r int = numbers.RandIntInRange(0,3)
						if r == 2 {
							(*nSeqPointer[p])[b] =  dna.ByteToBase(bases[3])
						} else {
							(*nSeqPointer[p])[b] =  dna.ByteToBase(bases[r])
						}
					case 'T':
						(*nSeqPointer[p])[b] =  dna.ByteToBase(bases[numbers.RandIntInRange(0,3)])
					}
					fmt.Println(nextFasta[p].Seq[b])
					fmt.Println(nextFasta[1].Seq[b])
					fmt.Println((*nSeqPointer[p])[b] == (*aSeqPointer[p])[b])
					fmt.Println(allFasta)
					fmt.Println((*nSeqPointer[p])[b])
				}
				
			}
		}
		/* 
		Have to make another for loop to update allFasta to be the same as nextFasta. 
		But it can't be updated inside the 3rd for loop. You want to get the entire population of 
		nextFasta completed before updating into allFasta, otherwise, some fasta at the end of
		nextFasta would be able to sample from the partially updated allFasta. 
		*/ 

		//fmt.Println(allFasta[0].Seq)
		//fmt.Println(nextFasta[0].Seq)

		//fmt.Println((*nPointer)[p].Seq )
		//fmt.Println(nextFasta[p].Seq)
		//for k := 0; k < popSize; k++ {
			// Update to allFasta
		//	(*aPointer)[k].Name = strconv.Itoa(k)
		//	*aSeqPointer[k] = *nSeqPointer[k]
		//}
		//fmt.Println(allFasta[p].Seq)

		

		//fmt.Println(allFasta[0].Seq)
	}

	return *aPointer
}