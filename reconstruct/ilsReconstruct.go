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
func IlsReconstructSeq(allPostProbs [][]float32, allRecons []pFasta.PFasta) pFasta.PFasta {
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
	var sum float32
	for pos := range allWeightedRecons[0].Seq {
		for idx := range weighted {
			sum += allWeightedRecons[idx].Seq[pos]
		}
		weightedSum.Seq[pos] = sum
	}

	return weightedSum

}