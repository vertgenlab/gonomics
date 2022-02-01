// Command Group: "Sequence Evolution & Reconstruction"

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

// likelihoodsToBaseUnbiased takes the un-normalized likelihoods for A, C, G, T as well
// and the probability
// threshold for when we will call it the mle base instead of N
// and gives back the reconstructed base for the hca
func likelihoodsToBaseUnbiased(likes []float64, probThreshold float64) dna.Base {
	var total, bestProb float64
	var i int
	var answer dna.Base = dna.N

	for i = range likes {
		total += likes[i]
	}

	for i = range likes {
		if likes[i]/total >= probThreshold && likes[i] > bestProb {
			bestProb = likes[i]
			answer = dna.Base(i)
		}
	}
	return answer
}

// likelihoodsToBaseHumanBias takes the un-normalized likelihoods for A, C, G, T as well
// as the index of the human base (0 for A, 1 for C, etc), and the probability
// threshold for when we will call it the mle base instead of the human base
// and gives back the reconstructed base for the hca
func likelihoodsToBaseHumanBias(likes []float64, humanBase dna.Base, probThreshold float64, nonHumanProbThreshold float64) dna.Base {
	var total, nonHumanTotal, bestProb float64
	var i int
	var answer dna.Base

	if humanBase == dna.Gap {
		answer = dna.N
	} else {
		answer = humanBase
	}

	for i = range likes {
		total += likes[i]
		if i != int(humanBase) {
			nonHumanTotal += likes[i]
		}
	}

	for i = range likes {
		if likes[i]/total >= probThreshold && likes[i] > bestProb && nonHumanTotal/total >= nonHumanProbThreshold {
			bestProb = likes[i]
			answer = dna.Base(i)
		}
	}
	return answer
}

func baseIsPresent(b dna.Base) bool {
	if dna.DefineBase(b) || b == dna.N {
		return true
	}
	return false
}

func hcaIsPresent(human, bonobo, chimp, gorilla, organutan dna.Base) bool {
	if baseIsPresent(human) && (baseIsPresent(bonobo) || baseIsPresent(chimp)) {
		return true
	}
	if (baseIsPresent(human) || baseIsPresent(bonobo) || baseIsPresent(chimp)) && (baseIsPresent(gorilla) || baseIsPresent(organutan)) {
		return true
	}
	return false
}

func reconHcaBase(root, humanNode, nodeToRecon *expandedTree.ETree, position int, probThreshold float64, nonHumanProbThreshold float64, humanBias bool) {
	var likelihoods []float64
	var nextBase dna.Base
	reconstruct.SetState(root, position)
	likelihoods = reconstruct.FixFc(root, nodeToRecon)
	if humanBias {
		nextBase = likelihoodsToBaseHumanBias(likelihoods, humanNode.Fasta.Seq[position], probThreshold, nonHumanProbThreshold)
	} else {
		nextBase = likelihoodsToBaseUnbiased(likelihoods, probThreshold)
	}
	nodeToRecon.Fasta.Seq = append(nodeToRecon.Fasta.Seq, nextBase)
}

func primateReconMle(inFastaFilename string, inTreeFilename string, humanBias bool, probThreshold float64, nonHumanProbThreshold float64, outputFastaFilename string) {
	var tree, humanNode, bonoboNode, chimpNode, gorillaNode, orangutanNode, hcaNode *expandedTree.ETree
	var err error
	var i int

	tree, err = expandedTree.ReadTree(inTreeFilename, inFastaFilename)
	exception.FatalOnErr(err)

	// roll call to make sure everyone is here and will need them later
	humanNode = expandedTree.FindNodeName(tree, "hg38")
	if humanNode == nil {
		log.Fatalf("Didn't find hg38 in the tree\n")
	}
	bonoboNode = expandedTree.FindNodeName(tree, "panPan2")
	if bonoboNode == nil {
		log.Fatalf("Didn't find panPan2 in the tree\n")
	}
	chimpNode = expandedTree.FindNodeName(tree, "panTro6")
	if chimpNode == nil {
		log.Fatalf("Didn't find panTro6 in the tree\n")
	}
	gorillaNode = expandedTree.FindNodeName(tree, "gorGor5")
	if gorillaNode == nil {
		log.Fatalf("Didn't find gorGor5 in the tree\n")
	}
	orangutanNode = expandedTree.FindNodeName(tree, "ponAbe3")
	if orangutanNode == nil {
		log.Fatalf("Didn't find ponAbe3 in the tree\n")
	}
	hcaNode = expandedTree.FindNodeName(tree, "hca")
	if hcaNode == nil {
		log.Fatalf("Didn't find hca in the tree\n")
	}

	for i = range humanNode.Fasta.Seq {
		if hcaIsPresent(humanNode.Fasta.Seq[i], bonoboNode.Fasta.Seq[i], chimpNode.Fasta.Seq[i], gorillaNode.Fasta.Seq[i], orangutanNode.Fasta.Seq[i]) {
			reconHcaBase(tree, humanNode, hcaNode, i, probThreshold, nonHumanProbThreshold, humanBias)
		} else {
			hcaNode.Fasta.Seq = append(hcaNode.Fasta.Seq, dna.Gap)
		}
	}
	fasta.Write(outputFastaFilename, []fasta.Fasta{*humanNode.Fasta, *hcaNode.Fasta})
}

func primateRecon(infile string, outfile string, messyToN bool) {
	records := fasta.Read(infile)
	output := append(records, reconstruct.PrimateRecon(records, messyToN))
	fasta.Write(outfile, output)
}

func usage() {
	fmt.Print(
		"primateRecon - Returns maximum likelihood sequence from an HBCGO primate alignment\n" +
			"Usage:\n" +
			"primateRecon input.fa output.fa\n" +
			"options:\n")
	flag.PrintDefaults()
}

func main() {
	var expectedNumArgs int = 2
	var messyToN *bool = flag.Bool("messyToN", false, "Sets messy bases to Ns in the output file.")
	var mleUnbiased *bool = flag.Bool("mleUnbiased", false, "Default is an N, unless a base passes the probThreshold threshold.")
	var mleHumanBiased *bool = flag.Bool("mleHumanBiased", false, "Default is the human base, unless the non-human bases, collectively and individually, pass nonHumanProbThreshold and probThreshold.")
	var tree *string = flag.String("mle", "", "Filename for newick tree with branch lengths.  Must have the anticipated assembly names and the hca.")
	var probThreshold *float64 = flag.Float64("probThreshold", 0.0, "The probability that a base other than human must pass to be considered a true change in the hca.")
	var nonHumanProbThreshold *float64 = flag.Float64("nonHumanProbThreshold", 0.0, "The sumation of all non-human bases must pass this threshold for a non-human base to be considered in the hca.")
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

	// some error check on flags provided
	if *messyToN && (*mleUnbiased || *mleHumanBiased) {
		log.Fatal("Error: -messyToN can not be used with mle estimates\n")
	}
	if (*tree == "") && (*mleUnbiased || *mleHumanBiased) {
		log.Fatal("Error: you need to provide a tree when using an mle estimate\n")
	}
	if *mleUnbiased && *mleHumanBiased {
		log.Fatal("Error: Can not do both a biased and unbiased mle estimate\n")
	}
	if (*probThreshold != 0 || *nonHumanProbThreshold != 0) && !(*mleUnbiased || *mleHumanBiased) {
		log.Fatal("Error: Can not use probability threshold flags without also using an mle estimate\n")
	}
	if *nonHumanProbThreshold != 0 && *mleUnbiased {
		log.Fatal("Error: Can not do a nonHumanProbThreshold when also doing an unbiased estimate\n")
	}

	if *mleUnbiased || *mleHumanBiased {
		// at this point we know that xor of the two mle flags is true, so we only pass one
		primateReconMle(inFile, *tree, *mleHumanBiased, *probThreshold, *nonHumanProbThreshold, outFile)
	} else {
		primateRecon(inFile, outFile, *messyToN)
	}
}
