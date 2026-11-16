// Package reconstruct has tools for reconstructing ancient genomes from the genomes of extant species using newick trees and fasta files.
package reconstruct

import (
	"log"

	"github.com/vertgenlab/gonomics/dna/pDna"
	"github.com/vertgenlab/gonomics/fasta"
	"github.com/vertgenlab/gonomics/fasta/pFasta"
	"github.com/vertgenlab/gonomics/wig"
)

// SumIlsSeqs takes in a list of posterior probabilities for n sequences and a list of n pFasta sequences
// and aggregates the pFasta sequences by weighing them by the corresponding posterior probability.
// Posterior probabilities should sum to 1 at each position across all tracks
func SumIlsSeqs(allPostProbs []map[string]wig.Wig, allRecons []pFasta.PFasta, precision float32) pFasta.PFasta {
	// checks
	if len(allPostProbs) != len(allRecons) {
		log.Fatalf("SumIlsSeqs requires the same number of postProbs (%v) as recons (%v)", len(allPostProbs), len(allRecons))
	}

	if len(allRecons) == 0 {
		log.Fatal("SumIlsSeqs requires at least one input sequence")
	}

	for idx, recon := range allRecons {
		currWig, ok := allPostProbs[idx][recon.Name]

		if !ok {
			log.Fatalf("No posterior probability WIG for %s", recon.Name)
		}

		if len(currWig.Values) != len(recon.Seq) {
			log.Fatalf("WIG and reconstruction lengths differ for %s", recon.Name)
		}

		if len(recon.Seq) != len(allRecons[0].Seq) {
			log.Fatal("Reconstructions do not have the same length")
		}
	}

	if len(allPostProbs[0][allRecons[0].Name].Values) != len(allRecons[0].Seq) {
		log.Fatalf("Reconstructions and posterior probabilies do not have the same length")
	}

	var lenRecon int

	lenRecon = len(allPostProbs[0][allRecons[0].Name].Values)
	if len(allRecons[0].Seq) < lenRecon {
		lenRecon = len(allRecons[0].Seq)
	}

	weightedSum := pFasta.PFasta{Name: "ilsRecon", Seq: make([]pDna.Float32Base, len(allRecons[0].Seq))}
	sum := pDna.Float32Base{A: 0, C: 0, G: 0, T: 0}
	zeroBase := pDna.Float32Base{A: 0, C: 0, G: 0, T: 0}

	var gapProbSum float32
	var nonGapProb float32

	for pos := range lenRecon {
		sum.A, sum.C, sum.G, sum.T = 0, 0, 0, 0

		// Sum weighted reconstructions at this position
		for idx, recon := range allRecons {
			prob := float32(allPostProbs[idx][recon.Name].Values[pos])
			sum = pDna.Sum(sum, pDna.Scale(recon.Seq[pos], prob))
		}

		if pDna.SumsToOne(sum, precision) {
			weightedSum.Seq[pos] = sum
			continue
		}

		gapProbSum = 0

		// Sum posterior probabilities of all gap-containing topologies
		for idx, recon := range allRecons {
			if pDna.IsGap(recon.Seq[pos]) {
				gapProbSum += float32(allPostProbs[idx][recon.Name].Values[pos])
			}
		}

		// Assign gap if combined posterior probability exceeds 0.5
		if gapProbSum > 0.5 {
			weightedSum.Seq[pos] = zeroBase
			continue
		}

		// Normalize by total probability of non-gap topologies
		nonGapProb = 1 - gapProbSum

		if nonGapProb <= 0 {
			weightedSum.Seq[pos] = zeroBase
			continue
		}

		weightedSum.Seq[pos] = pDna.Scale(sum, 1/nonGapProb)
		if !pDna.SumsToOne(weightedSum.Seq[pos], precision) {
			log.Fatalf("Normalized probabilities do not sum to one at position %d", pos)
		}

		if gapProbSum == 0 {
			log.Fatalf("Invalid nucleotide probability sum of zero at position %d", pos)
		}
	}

	return weightedSum
}

// SumIlsSeqsFasta converts topology-specific FASTA reconstructions into
// one-hot probabilistic sequences, then combines them using SumIlsSeqs.
func SumIlsSeqsFasta(allPostProbs []map[string]wig.Wig, allRecons []fasta.Fasta, precision float32) pFasta.PFasta {
	allRecons_pfa := make([]pFasta.PFasta, len(allRecons))

	for idx, record := range allRecons {
		allRecons_pfa[idx] = pFasta.FaToPfa(record, 0, len(record.Seq))
	}

	return SumIlsSeqs(allPostProbs, allRecons_pfa, precision)
}
