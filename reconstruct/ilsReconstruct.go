package reconstruct

import (
	"github.com/vertgenlab/gonomics/dna/pDna"
	"github.com/vertgenlab/gonomics/fasta"
	"github.com/vertgenlab/gonomics/fasta/pFasta"
	"github.com/vertgenlab/gonomics/wig"
	"log"
	"fmt"
	"strconv"
)

// IlsReconstructSeq takes in a list of posterior probabilities for n sequences and a list of n reconstructed pFasta sequences, and aggregates the pFasta sequences by weighing them by the corresponding posterior probability. Posterior probabilities should sum to 1 at each position.
func IlsReconstructSeq(allPostProbs []map[string]wig.Wig, allRecons []pFasta.PFasta, precision float32) pFasta.PFasta {
	// check that the posterior probabilities align with the given sequences
	if len(allPostProbs) != len(allRecons) {
		log.Fatalf("IlsReconstructSeq requires the same number of postProbs (%v) as recons (%v)", len(allPostProbs), len(allRecons))
	}

	//check that all the reconstructed sequences have the same length
	fmt.Println(len(allPostProbs))
	fmt.Println("allpostprobs 0")
	// for idx, item := range allPostProbs[0] {
	// 	fmt.Println("????")
	// }
	// fmt.Println(allPostProbs[0])
	fmt.Println("indexed into name")
	fmt.Println("allRecons name\t", len(allPostProbs[0][allRecons[0].Name].Values))
	fmt.Println("chr1\t", len(allPostProbs[0]["chr1"].Values))
	fmt.Println("v0\t", len(allPostProbs[0]["V0"].Values))

	w, ok := allPostProbs[0][allRecons[0].Name]
	fmt.Println("lookup ok:", ok, "name:", strconv.Quote(allRecons[0].Name), "len:", len(w.Values))
	if !ok {
		fmt.Println("available keys:")
		for k := range allPostProbs[0] {
			fmt.Println(strconv.Quote(k))
		}
	}

	for idx := 0; idx < len(allRecons)-1; idx++ {
		if len(allPostProbs[idx][allRecons[idx].Name].Values) != len(allPostProbs[idx+1][allRecons[idx+1].Name].Values) {
			log.Fatalf("Requested posterior probabilities do not have the same length.")
		}
		if len(allRecons[idx].Seq) != len(allRecons[idx+1].Seq) {
			log.Fatalf("Requested sequences do not have the same length.")
		}
	}

	// ans := pFasta.PFasta{Name: "chr1",
	// Seq: []pDna.Float32Base{
	// 	{
	// 		A: 0.23857,
	// 		C: 0.3323,
	// 		G: 0.44958,
	// 		T: 0.139448,
	// 	}}}
	
	// return ans
	
	// pre-allocate the allWeightedRecons slice
	var allWeightedRecons []pFasta.PFasta = make([]pFasta.PFasta, len(allRecons))

	// weight each sequence by the corresponding posterior probability
	for topologyIdx := range allPostProbs {
		// TODO: so theoretically, as I range through allPostProbs, it's not going to be
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
		}
		// } else {
		// 	log.Fatalf("This reconstruction returns a pDNA base that does not sum to 1 (%v) at %v", sum, pos)
		// }
	}

	return weightedSum

}

// IlsReconstructSeq takes in a list of posterior probabilities for n sequences and a list of n reconstructed pFasta sequences, and aggregates the pFasta sequences by weighing them by the corresponding posterior probability. Posterior probabilities should sum to 1 at each position.
func IlsReconstructSeqFasta(allPostProbs []map[string]wig.Wig, allRecons []fasta.Fasta, precision float32) pFasta.PFasta {
	allRecons_pfa := make([]pFasta.PFasta, 4)
	for idx, record := range allRecons {
		allRecons_pfa[idx] = pFasta.FaToPfa(record, 0, len(record.Seq))
	}

	return IlsReconstructSeq(allPostProbs, allRecons_pfa, precision)
}
