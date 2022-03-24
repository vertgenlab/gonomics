package simulate

import (
	"fmt"
	"math/rand"
	"strconv"

	"github.com/vertgenlab/gonomics/dna"
	"github.com/vertgenlab/gonomics/fasta"
	"github.com/vertgenlab/gonomics/numbers"
)

// SimulateAFS takes in a
func SimulateWrightFisher(popSize int, mutR float64, genTime int, genomeSize int, verbose int) []fasta.Fasta {
	curFasta := make([]fasta.Fasta, popSize)
	nextFasta := make([]fasta.Fasta, popSize)

	var r int
	var j, t, b, p, k int

	initialSeq := RandIntergenicSeq(0.5, genomeSize)

	// Copy the same fasta for the entire population
	for j = 0; j < popSize; j++ {
		curFasta[j].Name = fmt.Sprintf("Seq_%v", strconv.Itoa(j))
		curFasta[j].Seq = make([]dna.Base, genomeSize)
		copy(curFasta[j].Seq, initialSeq)

		nextFasta[j].Name = fmt.Sprintf("Seq_%v", strconv.Itoa(j))
		nextFasta[j].Seq = make([]dna.Base, genomeSize)
		copy(nextFasta[j].Seq, initialSeq)
	}

	// 1st loop through every generation
	for t = 1; t < genTime; t++ {

		// 2nd loop through every base
		for b = 0; b < genomeSize; b++ {

			// 3rd loop through every individual
			for p = 0; p < popSize; p++ {
				// Which chromosome from previous generation that pth individual inherits from?
				r = numbers.RandIntInRange(0, popSize)
				nextFasta[p].Seq[b] = curFasta[r].Seq[b] //randomize index of previous gen

				if verbose > 0 {
					fmt.Println(nextFasta)
				}

				// If random float generated is less than mutation rate (mutR), mutate the base
				if rand.Float64() < mutR {
					nextFasta[p].Seq[b] = changeBase(nextFasta[p].Seq[b])
				}
			}
		}
		/*
			Have to make another for loop to update allFasta to be the same as nextFasta.
			But it can't be updated inside the 3rd for loop. You want to get the entire population of
			nextFasta completed before updating into allFasta, otherwise, some fasta at the end of
			nextFasta would be able to sample from the partially updated allFasta.
		*/
		for k = 0; k < popSize; k++ {
			copy(curFasta, nextFasta)
		}
	}

	return curFasta
}
