package reconstruct

import (
	"log"

	"github.com/vertgenlab/gonomics/dna"
	"github.com/vertgenlab/gonomics/dna/pDna"
	"github.com/vertgenlab/gonomics/expandedTree"
	"github.com/vertgenlab/gonomics/fasta"
	"github.com/vertgenlab/gonomics/fasta/pFasta"
)

// ilsReconstructSeq takes in a list of posterior probabilities for n sequences and a list of n reconstructed pFasta sequences, and aggregates the pFasta sequences by weighing them by the corresponding posterior probability. Posterior probabilities should sum to 1 at each position.
func IlsReconstructSeq(allPostProbs [][]float32, allRecons []pFasta.PFasta, precision float32) pFasta.PFasta {
	if len(allPostProbs) != len(allRecons) {
		log.Fatalf("IlsReconstructSeq requires the same number of postProbs (%v) as recons (%v)", len(allPostProbs), len(allRecons))
	}
	for idx := range allRecons[:-1] {
		if len(allRecons[idx].Seq) != len(allRecons[idx+1].Seq) {
			log.Fatalf("Requested sequences do not have the same length.")
		}
	}
	
	var allWeightedRecons []pFasta.PFasta
	for topologyIdx := range postProbs {
		for pos, prob := range postProbs {
			allWeightedRecons[topologyIdx] = pDna.Scale(recons[topologyIdx].Seq[pos], prob)
		}
	}
	
	weightedSum := pFasta.PFasta{Name: recon1.Name, Seq = make([]pDna.Float32Base, len(recon1.Seq))}
	var sum pDna.Float32Base
	// iterate through each position of the sequence
	for pos := range allWeightedRecons[0].Seq {
		
		// sum over all the given recons at that position
		for _, recon := range allWeightedRecons {
			sum = pDna.Sum(sum, recon.Seq[pos])
		}
		
		if pDna.IsValid(sum, precision) {
			weightedSum.Seq[pos] = sum
		} else {
			log.Fatalf("This reconstruction returns an invalid pDNA base at %v", pos)
		}
		
	}

	return weightedSum

}