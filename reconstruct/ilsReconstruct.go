package reconstruct

import (
	"github.com/vertgenlab/gonomics/dna/pDna"
	"github.com/vertgenlab/gonomics/fasta/pFasta"
	"github.com/vertgenlab/gonomics/wig"
	"log"
)

// IlsReconstructSeq takes in a list of posterior probabilities for n sequences and a list of n reconstructed pFasta sequences, and aggregates the pFasta sequences by weighing them by the corresponding posterior probability. Posterior probabilities should sum to 1 at each position.
func IlsReconstructSeq(allPostProbs []map[string]wig.Wig, allRecons []pFasta.PFasta, precision float32) pFasta.PFasta {
	// check that the posterior probabilities align with the given sequences
	if len(allPostProbs) != len(allRecons) {
		log.Fatalf("IlsReconstructSeq requires the same number of postProbs (%v) as recons (%v)", len(allPostProbs), len(allRecons))
	}

	//check that all the reconstructed sequences have the same length
	for idx := 0; idx < len(allRecons)-1; idx++ {
		if len(allPostProbs[idx][allRecons[idx].Name].Values) != len(allPostProbs[idx+1][allRecons[idx+1].Name].Values) {
			log.Fatalf("Requested posterior probabilities do not have the same length.")
		}
		if len(allRecons[idx].Seq) != len(allRecons[idx+1].Seq) {
			log.Fatalf("Requested sequences do not have the same length.")
		}
	}

	// pre-allocate the allWeightedRecons slice
	var allWeightedRecons []pFasta.PFasta = make([]pFasta.PFasta, len(allRecons))

	// weight each sequence by the corresponding posterior probability
	for topologyIdx := range allPostProbs {
		allWeightedRecons[topologyIdx] = pFasta.PFasta{Name: allRecons[topologyIdx].Name, Seq: make([]pDna.Float32Base, len(allRecons[topologyIdx].Seq))}
		for pos, prob := range allPostProbs[topologyIdx][allRecons[topologyIdx].Name].Values {
			allWeightedRecons[topologyIdx].Seq[pos] = pDna.Scale(allRecons[topologyIdx].Seq[pos], float32(prob))
		}
	}

	// sum the weighted sequences
	weightedSum := pFasta.PFasta{Name: "ilsRecon", Seq: make([]pDna.Float32Base, len(allWeightedRecons[0].Seq))}
	sum := pDna.Float32Base{A: 0, C: 0, G: 0, T: 0}
	// iterate through each position of the sequence
	for pos := range allWeightedRecons[0].Seq {
		sum.A, sum.C, sum.G, sum.T = 0, 0, 0, 0
		// sum over all the given recons at that position
		for _, recon := range allWeightedRecons {
			sum = pDna.Sum(sum, recon.Seq[pos])
		}

		if pDna.SumsToOne(sum, precision) {
			weightedSum.Seq[pos] = sum
		} else {
			log.Fatalf("This reconstruction returns a pDNA base that does not sum to 1 at %v", pos)
		}
	}

	return weightedSum

}
