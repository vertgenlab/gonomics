package reconstruct

import (
	"github.com/vertgenlab/gonomics/dna"
	"github.com/vertgenlab/gonomics/expandedTree"
	"github.com/vertgenlab/gonomics/fasta"
	"github.com/vertgenlab/gonomics/genePred"
	"github.com/vertgenlab/gonomics/simulate"
	"log"
)

//returns the percentage accuracy by base returned by reconstruct of each node and of all reconstructed nodes combined.
//If option = 1 it will also return the percentage of first second and third bases of each codon that were correct
func ReconAccuracy(simFilename string, reconFilename string, leavesOnlyFile string, gpFilename string, option int) (accTotal map[string]float64, accBases map[string][]float64) {
	var accByBase map[string][]float64
	if option == 1 {
		accByBase = ReconAccuracyByBase(simFilename, reconFilename, gpFilename)
	}

	var allNodes string
	allNodes = "All Reconstructed Nodes"
	var found = false
	var leaf = false
	var exon = false
	var total float64
	total = 0.0
	var mistakes float64
	var leafMistakes float64
	var exonMistakes float64
	var nonCodingMistakes float64
	var exonBases float64
	var nonCodingBases float64
	sim := fasta.Read(simFilename)
	recon := fasta.Read(reconFilename)
	leaves := fasta.Read(leavesOnlyFile)
	genes := genePred.Read(gpFilename)

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
				//DEBUG: log.Printf("\n%s \n%s \n", dna.BasesToString(sim[i].Seq), dna.BasesToString(recon[j].Seq))
				for g := 0; g < len(genes); g++ {
					for k := 0; k < len(sim[i].Seq); k++ {
						exon, _ = simulate.CheckExon(genes[g], k)
						if exon == true {
							exonBases = exonBases + 1
						} else {
							nonCodingBases = nonCodingBases + 1
						}
						if sim[i].Seq[k] != recon[j].Seq[k] {
							if leaf == false {
								mistakes = mistakes + 1
							} else {
								leafMistakes = leafMistakes + 1
							}
							if exon == true {
								exonMistakes = exonMistakes + 1
							} else {
								nonCodingMistakes = nonCodingMistakes + 1
							}
						}
					}
				}
			}
		}
		if found == false {
			log.Fatal("Did not find all simulated sequences in reconstructed fasta.")
		}
		if leaf == false {
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

func ReconAccuracyByBase(simFilename string, reconFilename string, gpFilename string) map[string][]float64 {
	sim := fasta.Read(simFilename)
	recon := fasta.Read(reconFilename)
	genes := genePred.Read(gpFilename)
	var exonsProcessed map[int]bool
	var simCodons []simulate.CodonExt
	var reconCodons []simulate.CodonExt
	answer := make(map[string][]float64)
	var mistakes1 float64
	var mistakes2 float64
	var mistakes3 float64
	var total1 float64
	var total2 float64
	var total3 float64

	for s := 0; s < len(sim); s++ {
		for r := 0; r < len(recon); r++ {
			var percentage1 float64
			var percentage2 float64
			var percentage3 float64
			if sim[s].Name == recon[r].Name {
				for i := 0; i < len(sim[s].Seq); i++ {
					for g := 0; g < len(genes); g++ {
						inExon, exon := simulate.CheckExon(genes[g], i)
						if inExon {
							_, processed := exonsProcessed[exon]
							if !processed {
								simCodons = simulate.CreateCodons(simulate.BasesToBaseExt(sim[s].Seq), genes[g], exon)
								reconCodons = simulate.CreateCodons(simulate.BasesToBaseExt(recon[r].Seq), genes[g], exon)
								for c := 0; c < len(simCodons); c++ {
									if simCodons[c].Seq[0].SeqPos == reconCodons[c].Seq[0].SeqPos {
										for p := 0; p < 3; p++ {
											if p == 0 {
												total1 += 1
												if simCodons[c].Seq[p].Base != reconCodons[c].Seq[p].Base {
													mistakes1 += 1
												}
											} else if p == 1 {
												total2 += 1
												if simCodons[c].Seq[p].Base != reconCodons[c].Seq[p].Base {
													mistakes2 += 1
												}
											} else {
												total3 += 1
												if simCodons[c].Seq[p].Base != reconCodons[c].Seq[p].Base {
													mistakes3 += 1
												}
											}
										}
									}
								}
							}
						}
					}
				}
				percentage1 = (mistakes1 / total1) * 100
				percentage2 = (mistakes2 / total2) * 100
				percentage3 = (mistakes3 / total3) * 100
				answer[sim[s].Name] = append(answer[sim[s].Name], percentage1)
				answer[sim[s].Name] = append(answer[sim[s].Name], percentage2)
				answer[sim[s].Name] = append(answer[sim[s].Name], percentage3)
				log.Print(len(answer[sim[s].Name]))
			}
		}
	}
	return answer
}

//write assigned sequences at all nodes to a fasta file
func WriteTreeToFasta(tree *expandedTree.ETree, outFile string) {
	var fastas []fasta.Fasta
	nodes := expandedTree.GetTree(tree)

	for i := 0; i < len(nodes); i++ {
		fastas = append(fastas, *nodes[i].Fasta)
	}
	fasta.Write(outFile, fastas)
}

//write assigned sequences at leaf nodes to a fasta file
func WriteLeavesToFasta(tree *expandedTree.ETree, leafFile string) {
	var leafFastas []fasta.Fasta
	nodes := expandedTree.GetLeaves(tree)

	for i := 0; i < len(nodes); i++ {
		leafFastas = append(leafFastas, *nodes[i].Fasta)
	}
	fasta.Write(leafFile, leafFastas)
}

//calculate probability of switching from one base to another
func Prob(a int, b int, t float64) float64 {
	var p float64
	switch {
	case a > 3 || b > 3:
		p = 0
	case a == b:
		p = 1 - t
	default:
		p = t / 3
	}
	return p
}

//take in probability of all 4 bases return integer value of the most likely base
func Yhat(r []float64) int {
	var n float64
	n = 0
	var pos int
	for p, v := range r {
		if v > n {
			n = v
			pos = p
		}
	}
	return pos
}

func allZero(r []float64) bool {
	for _, v := range r {
		if v != 0 {
			return false
		}
	}
	return true
}

//set up Stored list for each node in the tree with probability of each base
func SetState(node *expandedTree.ETree, position int) {
	if node.Left != nil && node.Right != nil {
		SetState(node.Left, position)
		SetState(node.Right, position)
		if allZero(node.Left.Stored) && allZero(node.Right.Stored) {
			log.Fatal("no Stored values passed to internal node")
		}

		for i := 0; i < 4; i++ {
			sum := 0.0
			for j := 0; j < 4; j++ {
				for k := 0; k < 4; k++ {
					sum = sum + Prob(i, j, node.Left.BranchLength)*node.Left.Stored[j]*Prob(i, k, node.Right.BranchLength)*node.Right.Stored[k]
				}
			} //branch length (transition probability) times the probability of a base appearing given the context of the lower tree nodes
			node.Stored[i] = sum
		}
	} else if node.Left != nil {
		SetState(node.Left, position)
		if node.Left.Stored == nil {
			log.Fatal("no Stored values passed to internal node, left branch")
		}
		for i := 0; i < 4; i++ {
			sum := 0.0
			for j := 0; j < 4; j++ {
				sum = sum + Prob(i, j, node.Left.BranchLength)*node.Left.Stored[j]
			} //branch length (transition probability) times the probability of a base appearing given the context of the lower tree nodes
			node.Stored[i] = sum
		}
	} else if node.Right != nil {
		SetState(node.Right, position)
		if node.Right.Stored == nil {
			log.Fatal("no Stored values passed to internal node, right branch")
		}
		for i := 0; i < 4; i++ {
			sum := 0.0
			for k := 0; k < 4; k++ {
				sum = sum + Prob(i, k, node.Right.BranchLength)*node.Right.Stored[k]
			}
			node.Stored[i] = sum
		}
	} else if node.Right == nil && node.Left == nil {
		node.State = 4 // starts as N
		node.Stored = []float64{0, 0, 0, 0}

		if len(node.Fasta.Seq) <= position {
			log.Fatal("position specified is out of range of sequence \n")
		} else if len(node.Fasta.Seq) > position {
			node.State = int(node.Fasta.Seq[position])
			for i := 0; i < 4; i++ {
				if i == node.State {
					node.Stored[i] = 1
				} else {
					node.Stored[i] = 0
				}
			}
		}
	}
}

//Bubble up the tree using the memory of the previous nodes
func BubbleUp(node *expandedTree.ETree, prevNode *expandedTree.ETree, scrap []float64) {
	tot := 0.0
	scrapNew := []float64{0, 0, 0, 0}
	for i := 0; i < 4; i++ {
		sum := 0.0
		for j := 0; j < 4; j++ {
			for k := 0; k < 4; k++ {
				if prevNode.Up != nil {
					if prevNode == node.Left { //scrap is equal to one position of prevNode.Stored (Left or Right)
						sum = sum + Prob(i, j, node.Left.BranchLength)*Prob(i, k, node.Right.BranchLength)*scrap[j]*node.Right.Stored[k]
					} else if prevNode == node.Right {
						sum = sum + Prob(i, j, node.Left.BranchLength)*Prob(i, k, node.Right.BranchLength)*scrap[k]*node.Left.Stored[j]
					}
				} else if prevNode.Up == nil {
					sum = sum + Prob(i, j, node.Left.BranchLength)*Prob(i, k, node.Right.BranchLength)*node.Left.Stored[j]*node.Right.Stored[k]
				}
			}
		}
		scrapNew[i] = sum
	}
	if node.Up != nil {
		BubbleUp(node.Up, node, scrapNew)
	} else if node.Up == nil {
		tot = scrapNew[0] + scrapNew[1] + scrapNew[2] + scrapNew[3]
		node.Scrap = tot
	}
}

//fix each node and return the probabilities for each base at that site
func FixFc(root *expandedTree.ETree, node *expandedTree.ETree) []float64 {
	ans := []float64{0, 0, 0, 0}

	for i := 0; i < 4; i++ {
		scrap := []float64{0, 0, 0, 0} //checking one base at a time each time you call BubbleUp
		scrap[i] = node.Stored[i]
		if node.Up != nil {
			//Bubble up the tree using the memory of the previous node in relation to changing position taking in probabilities of bases
			//(node will be BubbleUp prevNode and node.Up will be the node being operated on)
			BubbleUp(node.Up, node, scrap) //node becomes PrevNode and scrap is set to one value of prevNode.Stored in BubbleUp
			ans[i] = root.Scrap            //root.Stored has previously assigned values (SetInternalState), you want to use whatever is returned by BubbleUp instead
		} else if node.Up == nil {
			ans[i] = root.Stored[i]
		}
	}

	return ans
}

//called by reconstructSeq.go on each base of the modern (leaf) seq. Loop over the nodes of the tree to return most probable base to the Fasta
func LoopNodes(root *expandedTree.ETree, position int) {
	internalNodes := expandedTree.GetBranch(root) //does this logic look like it will get leaves as well? We may want to use GetLeaves, which looks like it gets everything
	SetState(root, position)
	for k := 0; k < len(internalNodes); k++ {
		fix := FixFc(root, internalNodes[k])
		yHat := Yhat(fix)
		internalNodes[k].Fasta.Seq = append(internalNodes[k].Fasta.Seq, []dna.Base{dna.Base(yHat)}...)
	}
}
