// Package reconstruct has tools for reconstructing ancient genomes from the genomes of extant species using newick trees and fasta files.
package reconstruct

import (
	"log"

	"github.com/vertgenlab/gonomics/dna"
	"github.com/vertgenlab/gonomics/expandedTree"
	"github.com/vertgenlab/gonomics/fasta"
	"github.com/vertgenlab/gonomics/genePred"
	"github.com/vertgenlab/gonomics/simulate"
)

// ReconAccuracy returns the percentage accuracy by base returned by reconstruct of each node and of all reconstructed nodes combined.
// If calcBaseAcc = true it will also run ReconAccuracyByBase.
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

// WriteTreeToFasta writes assigned sequences at all nodes to a fasta file.
func WriteTreeToFasta(tree *expandedTree.ETree, outFile string) {
	var fastas []fasta.Fasta
	nodes := expandedTree.GetTree(tree)

	for i := 0; i < len(nodes); i++ {
		fastas = append(fastas, *nodes[i].Fasta)
	}
	fasta.Write(outFile, fastas)
}

// WriteLeavesToFasta writes assigned sequences at leaf nodes to a fasta file.
func WriteLeavesToFasta(tree *expandedTree.ETree, leafFile string) {
	var leafFastas []fasta.Fasta
	nodes := expandedTree.GetLeaves(tree)

	for i := 0; i < len(nodes); i++ {
		leafFastas = append(leafFastas, *nodes[i].Fasta)
	}
	fasta.Write(leafFile, leafFastas)
}

// mutationProbability calculate probability of switching from one base (a) to another (b) over branch length (t).
func mutationProbability(a int, b int, t float64) float64 {
	switch {
	case a > 3 || b > 3:
		return 0
	case a == b:
		return 1 - t
	default:
		return t / 3
	}
}

// LikelihoodsToBase returns the index of the most likely base for each position of a sequence. That position refers to a specific base.
func LikelihoodsToBase(likelihoods []float64, nonBiasBaseThreshold float64, biasBase dna.Base, highestProbThreshold float64) dna.Base {
	var highestProb, nonBiasBaseProb, total float64 = 0, 0, 0
	var answer dna.Base = biasBase
	for p, v := range likelihoods {
		total += v
		if dna.Base(p) != biasBase {
			nonBiasBaseProb += v
		}
		if v > highestProb {
			highestProb = v
			answer = dna.Base(p)
		}
	}
	if highestProb/total < highestProbThreshold {
		return dna.N
	}
	if nonBiasBaseProb/total < nonBiasBaseThreshold {
		return biasBase
	}
	return answer
}

// allZero checks if all values in a slice = 0 and returns a bool.
func allZero(r []float64) bool {
	for _, v := range r {
		if v != 0 {
			return false
		}
	}
	return true
}

// SetState calculates initial probabilities for all bases of the sequences for all nodes of the tree.
func SetState(node *expandedTree.ETree, position int) {
	var currNodeCounter, currLeftCounter, currRightCounter int
	var currSum float64

	if node.Left != nil && node.Right != nil {
		SetState(node.Left, position)
		SetState(node.Right, position)
		if allZero(node.Left.Stored) && allZero(node.Right.Stored) {
			log.Fatalf("Error: no Stored values passed to internal node.\n")
		}

		for currNodeCounter = range node.Stored {
			currSum = 0.0
			for currLeftCounter = range node.Left.Stored {
				for currRightCounter = range node.Right.Stored {
					currSum += mutationProbability(currNodeCounter, currLeftCounter, node.Left.BranchLength) * node.Left.Stored[currLeftCounter] * mutationProbability(currNodeCounter, currRightCounter, node.Right.BranchLength) * node.Right.Stored[currRightCounter]
				}
			} //branch length (transition probability) times the probability of a base appearing given the context of the lower tree nodes
			node.Stored[currNodeCounter] = currSum
		}
	} else if node.Left != nil {
		SetState(node.Left, position)
		if node.Left.Stored == nil {
			log.Fatalf("Error: no Stored values passed to internal node, left branch.\n")
		}
		for currNodeCounter = range node.Stored {
			currSum = 0.0
			for currLeftCounter = range node.Left.Stored {
				currSum += mutationProbability(currNodeCounter, currLeftCounter, node.Left.BranchLength) * node.Left.Stored[currLeftCounter]
			} //branch length (transition probability) times the probability of a base appearing given the context of the lower tree nodes
			node.Stored[currNodeCounter] = currSum
		}
	} else if node.Right != nil {
		SetState(node.Right, position)
		if node.Right.Stored == nil {
			log.Fatalf("Error: no Stored values passed to internal node, right branch.\n")
		}
		for currNodeCounter = range node.Stored {
			currSum = 0.0
			for currRightCounter = range node.Right.Stored {
				currSum += mutationProbability(currNodeCounter, currRightCounter, node.Right.BranchLength) * node.Right.Stored[currRightCounter]
			}
			node.Stored[currNodeCounter] = currSum
		}
	} else if node.Right == nil && node.Left == nil {
		node.State = 4 // starts as N
		node.Stored = []float64{0, 0, 0, 0}

		if len(node.Fasta.Seq) <= position {
			log.Fatal("position specified is out of range of sequence \n")
		} else {
			node.State = int(node.Fasta.Seq[position])
			for currNodeCounter = range node.Stored {
				if node.Fasta.Seq[position] == dna.N || node.Fasta.Seq[position] == dna.Gap {
					node.Stored[currNodeCounter] = 0.25
				} else if currNodeCounter == node.State {
					node.Stored[currNodeCounter] = 1
				} else {
					node.Stored[currNodeCounter] = 0
				}
			}
		}
	}
}

// bubbleUp calculates the final probabilities of all states in every position of the sequence at each internal node.
// using the stored values (initial probabilities from SetState) bubbleUp recursively calculates the
// probability on a child node based on both of the descendents of its ancestor. If a child node is a left child,
// bubbleUp uses the parent and right child's sequence information in Stored to compute a final probability
// of each of the base states at the left child, then passes those new probabilities up the tree to the root.
func bubbleUp(node *expandedTree.ETree, prevNode *expandedTree.ETree, scrap []float64) {
	var currNodeCounter, currLeftCounter, currRightCounter int
	var currSum, tot float64
	scrapNew := make([]float64, len(node.Stored))
	for currNodeCounter = range node.Stored {
		currSum = 0
		for currLeftCounter = range node.Left.Stored {
			for currRightCounter = range node.Right.Stored {
				if prevNode.Up != nil {
					if prevNode == node.Left { //scrap is equal to one position of prevNode.Stored (Left or Right)
						currSum += mutationProbability(currNodeCounter, currLeftCounter, node.Left.BranchLength) * mutationProbability(currNodeCounter, currRightCounter, node.Right.BranchLength) * scrap[currLeftCounter] * node.Right.Stored[currRightCounter]
					} else if prevNode == node.Right {
						currSum += mutationProbability(currNodeCounter, currLeftCounter, node.Left.BranchLength) * mutationProbability(currNodeCounter, currRightCounter, node.Right.BranchLength) * scrap[currRightCounter] * node.Left.Stored[currLeftCounter]
					}
				} else if prevNode.Up == nil {
					currSum += mutationProbability(currNodeCounter, currLeftCounter, node.Left.BranchLength) * mutationProbability(currNodeCounter, currRightCounter, node.Right.BranchLength) * node.Left.Stored[currLeftCounter] * node.Right.Stored[currRightCounter]
				}
			}
		}
		scrapNew[currNodeCounter] = currSum
	}
	if node.Up != nil {
		bubbleUp(node.Up, node, scrapNew)
	} else if node.Up == nil {
		tot = scrapNew[0] + scrapNew[1] + scrapNew[2] + scrapNew[3]
		node.Scrap = tot
	}
}

// FixFc passes a node to BubbleUp so that final base probabilities can be calculated from
// the initial values calculated in SetState.
func FixFc(root *expandedTree.ETree, node *expandedTree.ETree) []float64 {
	var currNodeCounter int
	ans := make([]float64, len(node.Stored))

	for currNodeCounter = range node.Stored {
		scrap := []float64{0, 0, 0, 0} //checking one base at a time each time you call BubbleUp
		scrap[currNodeCounter] = node.Stored[currNodeCounter]
		if node.Up != nil {
			//node will be bubbleUp prevNode and node.Up will be the node being operated on
			bubbleUp(node.Up, node, scrap)    //node becomes PrevNode and scrap is set to one value of prevNode.Stored in BubbleUp
			ans[currNodeCounter] = root.Scrap //root.Stored has previously assigned values (SetInternalState), you want to use whatever is returned by BubbleUp instead
		} else if node.Up == nil {
			ans[currNodeCounter] = root.Stored[currNodeCounter]
		}
	}

	return ans
}

// BaseExists determines, for a particular node at a particular position, whether a base exists (A, C, G, T, N)
// or not (node.Fasta.Seq[position] == dna.Gap). This information is stored as a bool in the field node.BasePresent.
// This function recursively searches descendent nodes and updates the BasePresent field.
func BaseExists(node *expandedTree.ETree, pos int) {
	var count int = 0
	if node.Left != nil && node.Right == nil {
		log.Fatalf("Error: tree is not a well-formed binary tree. For node: %v, left is not nil but right is nil.\n", node.Name)
	}
	if node.Right != nil && node.Left == nil {
		log.Fatalf("Error: tree is not a well-formed binary tree. For node: %v, right is not nil but left is nil.\n", node.Name)
	}
	if node.Up != nil && node.Up.BasePresent {
		count++
	}
	if node.Left != nil && node.Left.DescendentBasePresent {
		count++
	}
	if node.Right != nil && node.Right.DescendentBasePresent {
		count++
	}
	node.BasePresent = count >= 2
	if node.Left != nil {
		BaseExists(node.Left, pos)
	}
	if node.Right != nil {
		BaseExists(node.Right, pos)
	}
}

// DescendentBaseExists determines whether, for a given node, any descendent nodes contain a base.
// This bool output is stored in node.DescendentBasePresent.
// This function recursively updates this field for all descendent nodes of the input node.
func DescendentBaseExists(node *expandedTree.ETree, pos int) {
	if node.Left != nil && node.Right == nil {
		log.Fatalf("Error: tree is not a well-formed binary tree. For node: %v, left is not nil but right is nil.\n", node.Name)
	}
	if node.Right != nil && node.Left == nil {
		log.Fatalf("Error: tree is not a well-formed binary tree. For node: %v, right is not nil but left is nil.\n", node.Name)
	}
	if node.Left == nil && node.Right == nil { //if node is a leaf
		// a descendent base is present at a leaf if the sequence at this position is not a gap.
		node.DescendentBasePresent = node.Fasta.Seq[pos] != dna.Gap
	} else { //if we are not a leaf
		if node.Left != nil {
			DescendentBaseExists(node.Left, pos)
		}
		if node.Right != nil {
			DescendentBaseExists(node.Right, pos)
		}
		node.DescendentBasePresent = node.Left.DescendentBasePresent || node.Right.DescendentBasePresent
	}
}

// LoopNodes performs ancestral sequence reconstruction for all nodes of an input tree, specified by an input root node,
// at a user-specified alignment position.
// Options: the user may specify a 'biasLeafName'. When specified, the reconstruction of this sequence's immediate ancestor
// will be biased towards that descendent. The degree of this bias is controlled by the option 'nonBiasBaseThreshold'.
// The user may also specify a 'highestProbThreshold'. If the program is uncertain about ancestral reconstruction for a particular
// node, this option will allow LoopNodes to return an 'N' for that node instead.
func LoopNodes(root *expandedTree.ETree, position int, biasLeafName string, nonBiasBaseThreshold float64, highestProbThreshold float64) {
	var fix []float64
	var biasBase, answerBase dna.Base
	var biasParentName string
	var biasLeafNode *expandedTree.ETree

	if biasLeafName != "" {
		biasLeafNode = expandedTree.FindNodeName(root, biasLeafName)
		if biasLeafNode == nil {
			log.Fatalf("Didn't find %v in tree.\n", biasLeafName)
		}
		if biasLeafNode.Up == nil {
			log.Fatalf("Error: Bias reconstruction node was specified as the root node.")
		}
		biasParentName = biasLeafNode.Up.Name
	}

	internalNodes := expandedTree.GetBranch(root)
	SetState(root, position)
	DescendentBaseExists(root, position)
	BaseExists(root, position)
	for k := range internalNodes {
		fix = FixFc(root, internalNodes[k])

		if internalNodes[k].BasePresent {
			if biasParentName != "" && internalNodes[k].Name == biasParentName {
				biasBase = biasLeafNode.Fasta.Seq[position]
				answerBase = LikelihoodsToBase(fix, nonBiasBaseThreshold, biasBase, highestProbThreshold) //biased estimate
			} else {
				answerBase = LikelihoodsToBase(fix, 0, dna.N, highestProbThreshold) //unbiased estimate
			}
		} else {
			answerBase = dna.Gap
		}

		internalNodes[k].Fasta.Seq = append(internalNodes[k].Fasta.Seq, []dna.Base{answerBase}...)
	}
}
