package reconstruct

import (
	"github.com/vertgenlab/gonomics/dna/pDna"
	"github.com/vertgenlab/gonomics/fasta"
	"github.com/vertgenlab/gonomics/fasta/pFasta"
	"github.com/vertgenlab/gonomics/wig"
	"log"
	"strings"
)

// IlsReconstructSeq takes in a list of posterior probabilities for n sequences and a list of n reconstructed pFasta sequences, and aggregates the pFasta sequences by weighing them by the corresponding posterior probability. Posterior probabilities should sum to 1 at each position.
func IlsReconstructSeq(allPostProbs []map[string]wig.Wig, allRecons []pFasta.PFasta, precision float32) pFasta.PFasta {
	// checks
	if len(allPostProbs) != len(allRecons) {
		log.Fatalf("IlsReconstructSeq requires the same number of postProbs (%v) as recons (%v)", len(allPostProbs), len(allRecons))
	}

	for idx := 0; idx < len(allRecons)-1; idx++ {
		wigA, okA := allPostProbs[idx][allRecons[idx].Name]
		wigB, okB := allPostProbs[idx+1][allRecons[idx+1].Name]

		if !okA || !okB {
			keys := []string{}
			for k := range allPostProbs[idx] {
				keys = append(keys, k)
			}

			possible := strings.Join(keys, ", ")
			log.Fatalf("Reconstructions (%v, %v) do not have the same names as posterior probabilities (possible: %v)", okA, okB, possible)
		}

		if len(wigA.Values) != len(wigB.Values) {
			log.Fatalf("Posterior probabilities do not have the same length.")
		}

		if len(allRecons[idx].Seq) != len(allRecons[idx+1].Seq) {
			log.Fatalf("Reconstructions do not have the same length.")
		}
	}

	// weight each sequence by the corresponding posterior probability
	var allWeightedRecons []pFasta.PFasta = make([]pFasta.PFasta, len(allRecons))
	var topoName string
	var n int

	for topologyIdx := range allPostProbs {
		topoName = allRecons[topologyIdx].Name
		topoSeq := allRecons[topologyIdx].Seq

		allWeightedRecons[topologyIdx] = pFasta.PFasta{Name: topoName, Seq: make([]pDna.Float32Base, len(topoSeq))}
		topoProb := allPostProbs[topologyIdx][topoName].Values
		
		n = len(topoProb)
		if len(topoSeq) < n {
			n = len(topoSeq)
		}

		for pos := 0; pos < n; pos++ {
			allWeightedRecons[topologyIdx].Seq[pos] = pDna.Scale(topoSeq[pos], float32(topoProb[pos]))
		}
	}

	// Sum the weighted sequences
	// TODO: if prob = 0 at that position, default set it to V0
	// TODO: what if there is NOTHING at the position? 
	// should there be a boolean option to say if you want default {0.25, 0.25, 0.25, 0.25}
	
	n = len(allPostProbs[0][topoName].Values)
	if len(allRecons[0].Seq) < n {
		n = len(allRecons[0].Seq)
	}

	weightedSum := pFasta.PFasta{Name: "ilsRecon", Seq:  make([]pDna.Float32Base, len(allWeightedRecons[0].Seq)),}
	sum := pDna.Float32Base{A: 0, C: 0, G: 0, T: 0}
	zeroBase := pDna.Float32Base{A: 0, C: 0, G: 0, T: 0}

	for pos := range n {
		sum.A, sum.C, sum.G, sum.T = 0, 0, 0, 0
		// sum over all the given recons at that position
		for _, recon := range allWeightedRecons {
			sum = pDna.Sum(sum, recon.Seq[pos])
		}

		if pDna.SumsToOne(sum, precision) {
			weightedSum.Seq[pos] = sum
		} else
			for idx, recon := range allRecons {
				if pDna.IsGap(recon.Seq[pos]) {
					if allPostProbs[idx][allRecons[idx].Name].Values[pos] > 0.5 { // prob of gap greater than prob of non-gap
						weightedSum.Seq[pos] = zeroBase
					} else {
						weightedSum.Seq[pos] = pDna.Scale(sum , 1/(1-allPostProbs[idx][allRecons[idx].Name].Values[pos]))
					}
				}
			}
		}
	}
	
	return weightedSum
}

// IlsReconstructSeq takes in a list of posterior probabilities for n sequences and a list of n reconstructed pFasta sequences, and aggregates the pFasta sequences by weighing them by the corresponding posterior probability. Posterior probabilities should sum to 1 at each position.
func IlsReconstructSeqFasta(allPostProbs []map[string]wig.Wig, allRecons []fasta.Fasta, precision float32) pFasta.PFasta {
	allRecons_pfa := make([]pFasta.PFasta, len(allRecons))
	
	for idx, record := range allRecons {
		allRecons_pfa[idx] = pFasta.FaToPfa(record, 0, len(record.Seq))
	}
	
	return IlsReconstructSeq(allPostProbs, allRecons_pfa, precision)
}
