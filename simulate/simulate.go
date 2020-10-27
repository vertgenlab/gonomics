package simulate

import (
	"fmt"
	"github.com/vertgenlab/gonomics/dna"
	"github.com/vertgenlab/gonomics/expandedTree"
	"github.com/vertgenlab/gonomics/fasta"
	"github.com/vertgenlab/gonomics/genePred"
	"log"
	"math/rand"
)

var GC float64 = 0.42

// makes random gene with start and stop codon, must be length multiple of 3
func RandGene(name string, length int, GCcontent float64) []*fasta.Fasta {
	var AT float64
	AT = 1 - GCcontent
	seq := []dna.Base{dna.A, dna.T, dna.G}
	randLength := length - 6

	if length%3 != 0 {

		log.Fatal("length must be divisible by three")

	} else {

		for i := 0; i < randLength; i++ {
			r := rand.Float64()

			//cut-offs based on GC content of galGal6
			if r < GCcontent/2 {
				seq = append(seq, dna.G)
			} else if r < GCcontent {
				seq = append(seq, dna.C)
			} else if r < AT/2+GCcontent {
				seq = append(seq, dna.T)
			} else {
				seq = append(seq, dna.A)
			}

		}
	}

	r := rand.Float64()

	if r < 1/3 {
		seq = append(seq, dna.T, dna.A, dna.G)
	} else if r < 1/3*2 {
		seq = append(seq, dna.T, dna.G, dna.A)
	} else {
		seq = append(seq, dna.T, dna.A, dna.A)
	}

	record := fasta.Fasta{name, seq}
	var answer []*fasta.Fasta
	answer = append(answer, &record)
	return answer
}

//final function to run to simulate based off of the random gene and the tree
func Simulate(randSeqFilename string, root *expandedTree.ETree, gene string) {
	var rand1 []*fasta.Fasta

	rand1 = fasta.Read(randSeqFilename)
	root.Fasta = rand1[0]
	printSeqForNodes(root, rand1[0].Seq, gene)
}

// BLOSUM matrix for amino acid switching probabilities normalized to 0-1, unsure how it was calculated
var BLOSUM = [][]float64{[]float64{0.288590604, 0.03087248322, 0.03087248322, 0.02953020134, 0.02147651007, 0.0255033557, 0.04026845638, 0.07785234899, 0.01476510067, 0.04295302013, 0.05906040268, 0.04429530201, 0.01744966443, 0.02147651007, 0.02953020134, 0.08456375839, 0.04966442953, 0.005369127517, 0.01744966443, 0.06845637584, 0.0},
	[]float64{0.04457364341, 0.3449612403, 0.03875968992, 0.03100775194, 0.007751937984, 0.0484496124, 0.0523255814, 0.03294573643, 0.02325581395, 0.02325581395, 0.04651162791, 0.1201550388, 0.01550387597, 0.01744186047, 0.01937984496, 0.04457364341, 0.03488372093, 0.005813953488, 0.01744186047, 0.03100775194, 0.0},
	[]float64{0.05122494432, 0.04454342984, 0.3140311804, 0.08240534521, 0.008908685969, 0.03340757238, 0.04899777283, 0.06458797327, 0.03118040089, 0.02227171492, 0.03118040089, 0.05345211581, 0.01113585746, 0.01781737194, 0.02004454343, 0.06904231626, 0.04899777283, 0.004454342984, 0.01559020045, 0.02672605791, 0.0},
	[]float64{0.04104477612, 0.02985074627, 0.06902985075, 0.3973880597, 0.007462686567, 0.02985074627, 0.09141791045, 0.04664179104, 0.01865671642, 0.0223880597, 0.02798507463, 0.0447761194, 0.009328358209, 0.01492537313, 0.0223880597, 0.05223880597, 0.03544776119, 0.003731343284, 0.01119402985, 0.02425373134, 0.0},
	[]float64{0.06504065041, 0.0162601626, 0.0162601626, 0.0162601626, 0.4837398374, 0.01219512195, 0.0162601626, 0.0325203252, 0.008130081301, 0.04471544715, 0.06504065041, 0.02032520325, 0.0162601626, 0.02032520325, 0.0162601626, 0.0406504065, 0.03658536585, 0.00406504065, 0.01219512195, 0.05691056911, 0.0},
	[]float64{0.05588235294, 0.07352941176, 0.04411764706, 0.04705882353, 0.008823529412, 0.2147058824, 0.1029411765, 0.04117647059, 0.02941176471, 0.02647058824, 0.04705882353, 0.09117647059, 0.02058823529, 0.01470588235, 0.02352941176, 0.05588235294, 0.04117647059, 0.005882352941, 0.02058823529, 0.03529411765, 0.0},
	[]float64{0.05524861878, 0.04972375691, 0.04051565378, 0.09023941068, 0.007366482505, 0.06445672192, 0.2965009208, 0.0349907919, 0.02578268877, 0.02209944751, 0.03683241252, 0.07550644567, 0.01289134438, 0.01657458564, 0.02578268877, 0.05524861878, 0.03683241252, 0.005524861878, 0.01657458564, 0.03130755064, 0.0},
	[]float64{0.07827260459, 0.02294197031, 0.03913630229, 0.03373819163, 0.01079622132, 0.01889338731, 0.02564102564, 0.5101214575, 0.01349527665, 0.01889338731, 0.02834008097, 0.03373819163, 0.009446693657, 0.01619433198, 0.01889338731, 0.05128205128, 0.02968960864, 0.005398110661, 0.01079622132, 0.02429149798, 0.0},
	[]float64{0.04198473282, 0.04580152672, 0.0534351145, 0.03816793893, 0.007633587786, 0.03816793893, 0.0534351145, 0.03816793893, 0.3549618321, 0.02290076336, 0.03816793893, 0.04580152672, 0.01526717557, 0.03053435115, 0.01908396947, 0.04198473282, 0.02671755725, 0.007633587786, 0.0572519084, 0.02290076336, 0.0},
	[]float64{0.0471281296, 0.0176730486, 0.0147275405, 0.0176730486, 0.01620029455, 0.01325478645, 0.0176730486, 0.0206185567, 0.0088365243, 0.2709867452, 0.1678939617, 0.0235640648, 0.03681885125, 0.0441826215, 0.0147275405, 0.02503681885, 0.03976435935, 0.0058910162, 0.0206185567, 0.176730486, 0.0},
	[]float64{0.04453441296, 0.02429149798, 0.01417004049, 0.01518218623, 0.01619433198, 0.01619433198, 0.02024291498, 0.02125506073, 0.01012145749, 0.1153846154, 0.3755060729, 0.02530364372, 0.0495951417, 0.05465587045, 0.01417004049, 0.02429149798, 0.03340080972, 0.007085020243, 0.02226720648, 0.09615384615, 0.0},
	[]float64{0.05699481865, 0.1070811744, 0.0414507772, 0.0414507772, 0.008635578584, 0.05354058722, 0.07081174439, 0.04317789292, 0.0207253886, 0.02763385147, 0.04317789292, 0.2780656304, 0.01554404145, 0.01554404145, 0.02763385147, 0.05354058722, 0.03972366149, 0.00518134715, 0.01727115717, 0.03281519862, 0.0},
	[]float64{0.05220883534, 0.03212851406, 0.02008032129, 0.02008032129, 0.01606425703, 0.0281124498, 0.0281124498, 0.0281124498, 0.01606425703, 0.1004016064, 0.1967871486, 0.03614457831, 0.1606425703, 0.04819277108, 0.01606425703, 0.03614457831, 0.04016064257, 0.008032128514, 0.02409638554, 0.09236947791, 0.0},
	[]float64{0.03382663848, 0.01902748414, 0.01691331924, 0.01691331924, 0.01057082452, 0.01057082452, 0.01902748414, 0.02536997886, 0.01691331924, 0.06342494715, 0.1141649049, 0.01902748414, 0.02536997886, 0.3868921776, 0.01057082452, 0.02536997886, 0.02536997886, 0.01691331924, 0.088794926, 0.05496828753, 0.0},
	[]float64{0.05684754522, 0.02583979328, 0.02325581395, 0.03100775194, 0.01033591731, 0.02067183463, 0.03617571059, 0.03617571059, 0.01291989664, 0.02583979328, 0.03617571059, 0.04134366925, 0.01033591731, 0.01291989664, 0.4935400517, 0.04392764858, 0.03617571059, 0.002583979328, 0.01291989664, 0.03100775194, 0.0},
	[]float64{0.109947644, 0.04013961606, 0.05410122164, 0.04886561955, 0.01745200698, 0.03315881326, 0.05235602094, 0.06631762653, 0.01919720768, 0.02966841187, 0.04188481675, 0.05410122164, 0.01570680628, 0.02094240838, 0.02966841187, 0.219895288, 0.08202443281, 0.005235602094, 0.01745200698, 0.04188481675, 0.0},
	[]float64{0.07297830375, 0.03550295858, 0.04339250493, 0.03747534517, 0.01775147929, 0.02761341223, 0.03944773176, 0.04339250493, 0.01380670611, 0.05325443787, 0.0650887574, 0.04536489152, 0.01972386588, 0.02366863905, 0.02761341223, 0.09270216963, 0.2465483235, 0.005917159763, 0.01775147929, 0.07100591716, 0.0},
	[]float64{0.0303030303, 0.02272727273, 0.01515151515, 0.01515151515, 0.007575757576, 0.01515151515, 0.02272727273, 0.0303030303, 0.01515151515, 0.0303030303, 0.05303030303, 0.02272727273, 0.01515151515, 0.06060606061, 0.007575757576, 0.02272727273, 0.02272727273, 0.4924242424, 0.06818181818, 0.0303030303, 0.0},
	[]float64{0.04049844237, 0.02803738318, 0.02180685358, 0.01869158879, 0.009345794393, 0.02180685358, 0.02803738318, 0.02492211838, 0.04672897196, 0.04361370717, 0.06853582555, 0.03115264798, 0.01869158879, 0.1308411215, 0.01557632399, 0.03115264798, 0.02803738318, 0.02803738318, 0.3177570093, 0.04672897196, 0.0},
	[]float64{0.06995884774, 0.0219478738, 0.01646090535, 0.01783264746, 0.01920438957, 0.01646090535, 0.02331961591, 0.02469135802, 0.008230452675, 0.1646090535, 0.1303155007, 0.02606310014, 0.03155006859, 0.03566529492, 0.01646090535, 0.0329218107, 0.04938271605, 0.00548696845, 0.02057613169, 0.268861454, 0.0},
	[]float64{0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0}}

//choose base based off of random float, takes in GC content
func chooseRandomBase(GCcontent float64) dna.Base {
	var base dna.Base
	var AT float64
	AT = 1 - GCcontent

	r := rand.Float64()

	if r < GCcontent/2 {
		base = dna.G
	} else if r < GCcontent {
		base = dna.C
	} else if r < AT/2+GCcontent {
		base = dna.T
	} else {
		base = dna.A
	}

	return base
}

//calls chooseRandomBase and loops until a different base than the original is found
func changeBase(originalBase dna.Base) dna.Base {
	newBase := chooseRandomBase(GC)

	for newBase == originalBase {
		newBase = chooseRandomBase(GC)
	}
	return newBase
}

//mutate base given random float, whether it's mutated is dependent on branchLength
func mutateBase(b dna.Base, branchLength float64) dna.Base {
	r := rand.Float64()

	var base dna.Base

	if branchLength == 0 {
		base = b
	} else if r < branchLength {
		base = changeBase(b)
	} else {
		base = b
	}

	return base
}

//mutate sequence taking BLOSUM probabilites and gene structure into account
func MutateSeq(inputSeq []dna.Base, branchLength float64, gene string) []dna.Base {
	var originalBase dna.Base
	var newBase dna.Base
	var originalCodons []*dna.Codon
	var newCodons []*dna.Codon
	var originalAmAc dna.AminoAcid
	var newAmAc dna.AminoAcid
	var newSequence []dna.Base
	var basesProcessed int
	var geneRecord []*genePred.GenePred

	seq := copySeq(inputSeq)
	geneRecord = genePred.Read(gene)

	//p will be inaccurate if in a coding sequence. basesProcessed variable will reflect how many bases have gone through codon simulation.
	for g := 0; g < len(geneRecord); g++ {
		for p := 0; p < len(seq); p++ {
			overlapExon, thisExon := CheckExon(geneRecord[g], p)
			log.Printf("position: %v", p)

			if overlapExon == false {
				newBase = mutateBase(seq[p], branchLength)
				newSequence = append(newSequence, newBase)
				//DEBUG:fmt.Printf("newSequence to position %v: %s\n", p+1, dna.BasesToString(newSequence))
			} else {
				if (geneRecord[g].ExonEnds[thisExon]-geneRecord[g].ExonStarts[thisExon]+1)%3 != 0 {
					log.Fatal("sequence length must be divisible by three")
				} else {
					codonNum := (geneRecord[g].ExonEnds[thisExon] - geneRecord[g].ExonStarts[thisExon] + 1) / 3
					//fmt.Printf("codonNum: %v\n position: %v\n", codonNum, p+1)

					for i := 0; i < codonNum; i++ {
						originalCodons = dna.BasesToCodons(seq)

						for j := 0; j < 3; j++ {
							originalBase = originalCodons[i].Seq[j]
							var thisCodon []dna.Base

							if i == 0 && thisExon == 0 {
								newBase = originalBase //cannot change start codon
							} else if i == codonNum-1 && thisExon == len(geneRecord[g].ExonStarts) { //if we are on the last codon of the last exon
								r := rand.Float64()
								if j == 0 { //first position is only ever a T
									originalCodons[i].Seq[j] = dna.T
								} else if j == 1 { //second position can either be an A or G
									if r < 0.66 {
										originalCodons[i].Seq[j] = dna.A
									} else {
										originalCodons[i].Seq[j] = dna.G
									}
								} else if j == 2 { //last position can either be A or G, but if previous position is G it cannot be G again
									if originalCodons[i].Seq[j-1] == dna.G {
										originalCodons[i].Seq[j] = dna.A
									} else {
										if r < 0.5 {
											originalCodons[i].Seq[j] = dna.A
										} else {
											originalCodons[i].Seq[j] = dna.G
										}
									}
								}
							} else {
								newBase = mutateBase(originalBase, branchLength)

								if j == 0 {
									thisCodon = append(thisCodon, newBase)
									thisCodon = append(thisCodon, originalCodons[i].Seq[j+1])
									thisCodon = append(thisCodon, originalCodons[i].Seq[j+2])
								} else if j == 1 {
									thisCodon = append(thisCodon, originalCodons[i].Seq[j-1])
									thisCodon = append(thisCodon, newBase)
									thisCodon = append(thisCodon, originalCodons[i].Seq[j+1])
								} else {
									thisCodon = append(thisCodon, originalCodons[i].Seq[j-2])
									thisCodon = append(thisCodon, originalCodons[i].Seq[j-1])
									thisCodon = append(thisCodon, newBase)
								}

								newCodons = dna.BasesToCodons(thisCodon)
								originalAmAc = dna.TranslateCodon(originalCodons[i])
								newAmAc = dna.TranslateCodon(newCodons[0])

								prob := BLOSUM[originalAmAc][newAmAc]
								r := rand.Float64()

								if r < prob {
									originalCodons[i].Seq[j] = newBase
								} else {
									originalCodons[i].Seq[j] = originalBase
								}
							}
							newSequence = append(newSequence, originalCodons[i].Seq[j])
							//DEBUG:fmt.Printf("newSequence @%v: %s\n", p+1, dna.BasesToString(newSequence))
							basesProcessed++
							fmt.Printf("basesProcessed: %v\n p: %v\n", basesProcessed, p)
						}
					}
					p += 2 + (3 * (codonNum - 1)) //prevents looping through already processed bases of the codon which are handled within j loop
				}
			}
		}
	}

	return newSequence
}

//func getOverlapCDS(gtf map[string]*gtf.Gene, pos int) (cds *gtf.CDS, start bool) { //only works on pos strand genes
//	for _, gene := range gtf {
//		for _, transcript := range gene.Transcripts {
//			for exonNum, exon := range transcript.Exons {
//				if exon.Cds != nil {
//					if exon.Cds.Start-1 <= pos && pos <= exon.Cds.End-1 {
//						return exon.Cds, exonNum == 0
//					}
//				}
//			}
//		}
//	}
//	return nil, false
//}
//
//func getOverlap5UTR(gtf map[string]*gtf.Gene, pos int) (UTR *gtf.FiveUTR, start bool) {
//	for _, gene := range gtf {
//		for _, transcript := range gene.Transcripts {
//			for exonNum, exon := range transcript.Exons {
//				if exon.FiveUtr != nil {
//					if exon.FiveUtr.Start-1 <= pos && pos <= exon.FiveUtr.End-1 {
//						return exon.FiveUtr, exonNum == 0
//					}
//				}
//			}
//		}
//	}
//	return nil, false
//}
//
//func getOverlap3UTR(gtf map[string]*gtf.Gene, pos int) (UTR *gtf.ThreeUTR, start bool) {
//	for _, gene := range gtf {
//		for _, transcript := range gene.Transcripts {
//			for exonNum, exon := range transcript.Exons {
//				if exon.ThreeUtr != nil {
//					if exon.ThreeUtr.Start-1 <= pos && pos <= exon.ThreeUtr.End-1 {
//						return exon.ThreeUtr, exonNum == 0
//					}
//				}
//			}
//		}
//	}
//	return nil, false
//}
//

func CheckExon(gene *genePred.GenePred, position int) (bool, int) {
	var answer bool
	var exonNum int
	for i := 0; i < len(gene.ExonStarts); i++ {
		if position >= gene.ExonStarts[i] && position <= gene.ExonEnds[i] {
			answer = true
			exonNum = i
		} else {
			answer = false
			exonNum = i
		}
	}
	return answer, exonNum
}

//make a slice and a copy of that list of an original sequence so the sequence can be assigned to a node and then mutated
func copySeq(seq []dna.Base) []dna.Base {
	original := make([]dna.Base, len(seq))
	copy(original, seq)
	return original
}

//make fastas based off of node and random sequence
func printSeqForNodes(node *expandedTree.ETree, sequence []dna.Base, gene string) {
	var length float64
	var seq []dna.Base
	var seqFasta fasta.Fasta
	//var fastaFinal []*fasta.Fasta

	length = node.BranchLength
	seq = MutateSeq(sequence, length, gene)

	seqFasta = fasta.Fasta{node.Name, seq}
	node.Fasta = &seqFasta
	//fastaFinal = append(fastaFinal, &seqFasta)
	if node.Left != nil && node.Right != nil {
		printSeqForNodes(node.Right, seq, gene)
		//fastaFinal = append(fastaFinal, b...)
		printSeqForNodes(node.Left, seq, gene)
		//fastaFinal = append(fastaFinal, a...)
	}
}

func RemoveAncestors(filename string, tree *expandedTree.ETree, outputFilename string) {
	var fastas []*fasta.Fasta
	var newFastas []*fasta.Fasta
	var outFile string

	fastas = fasta.Read(filename)

	leaves := expandedTree.GetLeaves(tree)
	for i := 0; i < len(fastas); i++ {
		for j := 0; j < len(leaves); j++ {
			if fastas[i].Name == leaves[j].Name {
				newFastas = append(newFastas, fastas[i])
			}
		}
	}
	outFile = outputFilename
	fasta.Write(outFile, newFastas)
}
