package simulate

import (
	"fmt"
	"github.com/vertgenlab/gonomics/fasta"
	"github.com/vertgenlab/gonomics/dna"
	"github.com/vertgenlab/gonomics/numbers"
)

// GoSimulateAFS takes in a
func GoSimulateAFS (popSize uint, mutRate int, genTime uint, genomeSize int) []fasta.Fasta {
	var initialGenome fasta.Fasta
	initialSeq := make([]dna.Base, genomeSize)
	var bases = []byte{'A', 'C', 'G', 'T'}

	// Randomly generate the initial genome sequence
	for i := 0; i < genomeSize; i++ {
		initialSeq[i] = dna.ByteToBase(bases[numbers.RandIntInRange(0,4)])
	}
	

	initialGenome.Name = "Adam"
	initialGenome.Seq = initialSeq

	fmt.Println(initialGenome.Name, initialGenome.Seq)
	answer := []fasta.Fasta{initialGenome}
	return answer
}