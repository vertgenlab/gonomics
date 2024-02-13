package reconstruct

import (
	"github.com/vertgenlab/gonomics/fasta"
	"github.com/vertgenlab/gonomics/genePred"
	"github.com/vertgenlab/gonomics/simulate"
	"log"
)

// ReconAccuracy returns the percentage accuracy by base returned by reconstruct of each node and of all reconstructed nodes combined.
// If calcBaseAcc = true it will also run ReconAccuracyByBase.
// accTotal return is a map of node names to accuracy float.
// accBases returns a map of node names to a []float64, corresponding to the accuracy for each base, A, C, G, T.
func ReconAccuracy(simFilename string, reconFilename string, leavesOnlyFile string, gpFilename string, calcBaseAcc bool) (accTotal map[string]float64, accBases map[string][]float64) {
	var accByBase map[string][]float64
	if calcBaseAcc {
		accByBase = ReconAccuracyByBase(simFilename, reconFilename, gpFilename)
	}

	var allNodes string
	allNodes = "All Reconstructed Nodes"
	var found = false
	var leaf = false
	var exon = false
	var total = 0.0
	var mistakes, leafMistakes, exonMistakes, nonCodingMistakes float64
	var exonBases, nonCodingBases float64
	sim := fasta.Read(simFilename)
	recon := fasta.Read(reconFilename)
	leaves := fasta.Read(leavesOnlyFile)
	var genes = make([]genePred.GenePred, 0)
	if gpFilename != "" {
		genes = genePred.Read(gpFilename)
	}

	answer := make(map[string]float64)

	for i := 0; i < len(sim); i++ {
		exonBases = 0.0
		nonCodingBases = 0.0
		mistakes = 0.0
		exonMistakes = 0.0
		nonCodingMistakes = 0.0
		found = false
		for j := 0; j < len(recon); j++ {
			if sim[i].Name == recon[j].Name {
				leaf = false
				for l := 0; l < len(leaves); l++ {
					if recon[j].Name == leaves[l].Name {
						leaf = true
					}
				}
				found = true
				for g := 0; g < len(genes); g++ {
					for k := 0; k < len(sim[i].Seq); k++ {
						exon, _ = simulate.CheckExon(genes[g], k)
						if exon {
							exonBases = exonBases + 1
						} else {
							nonCodingBases = nonCodingBases + 1
						}
						if sim[i].Seq[k] != recon[j].Seq[k] {
							if !leaf {
								mistakes = mistakes + 1
							} else {
								leafMistakes = leafMistakes + 1
							}
							if exon {
								exonMistakes = exonMistakes + 1
							} else {
								nonCodingMistakes = nonCodingMistakes + 1
							}
						}
					}
				}
			}
		}
		if !found {
			log.Fatal("Did not find all simulated sequences in reconstructed fasta.")
		}
		if !leaf {
			accuracy := mistakes / float64(len(sim[i].Seq)) * 100.0
			//DEBUG: fmt.Printf("tot: %f, len(sim): %f, len(sim[0].Seq): %f \n", total, float64(len(sim)), float64(len(sim[0].Seq)))
			acc := 100 - accuracy
			answer[sim[i].Name] = acc
			total = total + mistakes
		} else {
			leafAccuracy := leafMistakes / float64(len(sim[i].Seq)) * 100.0
			lAcc := 100 - leafAccuracy
			answer[sim[i].Name+"(leaf)"] = lAcc
		}
		exonAccuracy := exonMistakes / exonBases * 100.0
		nonCodingAccuracy := nonCodingMistakes / nonCodingBases * 100.0
		//DEBUG: log.Printf("Node: %s, nonCodingMistakes: %f, nonCodingBases: %f, exonMistakes: %f, exonBases: %f", sim[i].Name, nonCodingMistakes, nonCodingBases, exonMistakes, exonBases)
		eAcc := 100 - exonAccuracy
		nCAcc := 100 - nonCodingAccuracy
		answer[sim[i].Name+" exon"] = eAcc
		answer[sim[i].Name+" nonCoding"] = nCAcc
	}
	accuracy := total / (float64(len(sim)-len(leaves)) * float64(len(sim[0].Seq))) * 100.0
	//DEBUG: fmt.Printf("tot: %f, acc: %f, len(sim): %f, len(leaves): %f, len(sim[0].Seq): %f \n", total, accuracy, float64(len(sim)), float64(len(leaves)), float64(len(sim[0].Seq)))
	acc := 100 - accuracy
	answer[allNodes] = acc

	return answer, accByBase
}

// ReconAccuracyByBase will run if the calcBaseAcc argument of ReconAccuracy = true.
// This function calculates the percentage of first, second and third bases in codons that were correct.
func ReconAccuracyByBase(simFilename string, reconFilename string, gpFilename string) map[string][]float64 {
	sim := fasta.Read(simFilename)
	rec := fasta.Read(reconFilename)
	recon := fasta.ToMap(rec)
	genes := genePred.Read(gpFilename)
	answer := make(map[string][]float64)

	for s := 0; s < len(sim); s++ {
		var percentage0, percentage1, percentage2 float64
		var mistakes0, mistakes1, mistakes2 float64
		var total0, total1, total2 float64
		rSeq, ok := recon[sim[s].Name]
		if ok {
			for i := 0; i < len(sim[s].Seq); i++ {
				for g := 0; g < len(genes); g++ {
					inExon, exon := simulate.CheckExon(genes[g], i)
					if inExon {
						loc := calcLocationInCodon(genes[g], exon, i)
						switch loc {
						case 0:
							total0 += 1
							if sim[s].Seq[i] != rSeq[i] {
								mistakes0 += 1
							}
						case 1:
							total1 += 1
							if sim[s].Seq[i] != rSeq[i] {
								mistakes1 += 1
							}
						case 2:
							total2 += 1
							if sim[s].Seq[i] != rSeq[i] {
								mistakes2 += 1
							}
						}
					}
				}
			}
			percentage0 = (mistakes0 / total0) * 100
			percentage1 = (mistakes1 / total1) * 100
			percentage2 = (mistakes2 / total2) * 100
			answer[sim[s].Name] = append(answer[sim[s].Name], 100-percentage0)
			answer[sim[s].Name] = append(answer[sim[s].Name], 100-percentage1)
			answer[sim[s].Name] = append(answer[sim[s].Name], 100-percentage2)
		} else {
			log.Panicf("Cannot find a reconstructed sequence match for simulated sequence: %s.", sim[s].Name)
		}
	}
	return answer
}

// calcLocationInCodon returns what position in a codon a given base inhabits (0, 1, 2).
func calcLocationInCodon(gene genePred.GenePred, exon int, position int) int {
	var positionInCodon int

	location := position - gene.ExonStarts[exon] + gene.ExonFrames[exon]
	positionInCodon = location % 3
	return positionInCodon
}
