package main

import (
	"flag"
	"fmt"
	"github.com/vertgenlab/gonomics/dna"
	"github.com/vertgenlab/gonomics/exception"
	"github.com/vertgenlab/gonomics/expandedTree"
	"github.com/vertgenlab/gonomics/fasta"
	"github.com/vertgenlab/gonomics/reconstruct"
	"log"
)

func mouseReconMraMle(inFile string, outFile string, treeFileName string, probThreshold float64, nonBiasProbThreshold float64) {
	var tree, mouseNode, ratNode, hamsterNode, squirrelNode, mraNode *expandedTree.ETree
	var err error
	var i int

	tree, err = expandedTree.ReadTree(treeFileName, inFile)
	exception.FatalOnErr(err)

	//roll call
	mouseNode = expandedTree.FindNodeName(tree, "mm10")
	if mouseNode == nil {
		log.Fatalf("Didn't find mm10 in the tree\n")
	}
	ratNode = expandedTree.FindNodeName(tree, "rn7")
	if ratNode == nil {
		log.Fatalf("Didn't find rn7 in the tree\n")
	}
	hamsterNode = expandedTree.FindNodeName(tree, "criGriChoV2")
	if hamsterNode == nil {
		log.Fatalf("Didn't find criGriChoV2 in the tree\n")
	}
	squirrelNode = expandedTree.FindNodeName(tree, "speTri2")
	if squirrelNode == nil {
		log.Fatalf("Didn't find speTri2 in the tree\n")
	}
	mraNode = expandedTree.FindNodeName(tree, "mra")
	if mraNode == nil {
		log.Fatalf("Didn't find mra in the tree\n")
	}

	for i = range mouseNode.Fasta.Seq {
		if mraIsPresent(mouseNode.Fasta.Seq[i], ratNode.Fasta.Seq[i], hamsterNode.Fasta.Seq[i], squirrelNode.Fasta.Seq[i]) {
			reconMraBase(tree, mouseNode, mraNode, i, probThreshold, nonBiasProbThreshold)
		} else {
			mraNode.Fasta.Seq = append(mraNode.Fasta.Seq, dna.Gap)
		}
	}
	fasta.Write(outFile, []fasta.Fasta{*mouseNode.Fasta, *ratNode.Fasta, *hamsterNode.Fasta, *squirrelNode.Fasta, *mraNode.Fasta})

}

func mraIsPresent(mouse, rat, hamster, squirrel dna.Base) bool {
	if baseIsPresent(mouse) && baseIsPresent(rat) {
		return true
	}
	if (baseIsPresent(mouse) || baseIsPresent(rat)) && (baseIsPresent(hamster) || baseIsPresent(squirrel)) {
		return true
	}
	return false
}

func baseIsPresent(b dna.Base) bool {
	if dna.DefineBase(b) || b == dna.N {
		return true
	}
	return false
}

func reconMraBase(root, mouseNode, nodeToRecon *expandedTree.ETree, position int, probThreshold float64, nonBiasProbThreshold float64) {
	var likelihoods []float64
	var nextBase dna.Base
	reconstruct.SetState(root, position)
	likelihoods = reconstruct.FixFc(root, nodeToRecon)
	nextBase = likelihoodToBaseBias(likelihoods, mouseNode.Fasta.Seq[position], probThreshold, nonBiasProbThreshold)
	nodeToRecon.Fasta.Seq = append(nodeToRecon.Fasta.Seq, nextBase)
}

func likelihoodToBaseBias(likes []float64, biasBase dna.Base, probThreshold float64, nonBiasProbThreshold float64) dna.Base {
	var total, nonBiasTotal, bestProb float64
	var i int
	var answer dna.Base

	if biasBase == dna.Gap {
		answer = dna.N
	} else {
		answer = biasBase
	}

	for i = range likes {
		total += likes[i]
		if i != int(biasBase) {
			nonBiasTotal += likes[i]
		}
	}

	for i = range likes {
		if likes[i]/total >= probThreshold && likes[i] > bestProb && nonBiasTotal/total >= nonBiasProbThreshold {
			bestProb = likes[i]
			answer = dna.Base(i)
		}
	}
	return answer
}

func usage() {
	fmt.Print(
		"mouseRecon - Returns maximum likelihood ancestral sequences from a Mouse-Rat-ChineseHamster-Squirrel multiple alignment in multiFa format.\n" +
			"Usage:\n" +
			"mouseRecon input.fa output.fa\n" +
			"options:\n")
	flag.PrintDefaults()
}

func main() {
	var expectedNumArgs int = 2
	var tree *string = flag.String("mleTree", "", "Filename for newick tree with branch lengths. Must have the anticipated assembly names and the desired ancestral node labeled.")
	var probThreshold *float64 = flag.Float64("probThreshold", 0.0, "The minimum probability assigned to a base required for it to be called.")
	var nonBiasProbThreshold *float64 = flag.Float64("nonBiasProbThreshold", 0.0, "The summation of all base probabilities (other than the biased species) must pass this threshold for a non-biased species base to be considered in the ancestral state.")
	//var mleMraUnbiased *bool = flag.Bool("mleMraUnbiased", false, "Calculates the unbiased mouse-rat ancestor estimation. Currently unsupported (future feature).")
	var mleMraMouseBiased *bool = flag.Bool("mleMraMouseBiased", true, "Calculates the mouse-biased mouse-rat ancestor estimation (default behavior, only behavior currently supported).")

	flag.Usage = usage
	log.SetFlags(log.Ldate | log.Ltime | log.Lshortfile)
	flag.Parse()

	if len(flag.Args()) != expectedNumArgs {
		flag.Usage()
		log.Fatalf("Error: expecting %d arguments, but got %d\n",
			expectedNumArgs, len(flag.Args()))
	}

	inFile := flag.Arg(0)
	outFile := flag.Arg(1)

	if *mleMraMouseBiased {
		mouseReconMraMle(inFile, outFile, *tree, *probThreshold, *nonBiasProbThreshold)
	}
}
