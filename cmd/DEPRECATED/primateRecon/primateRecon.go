// Command Group: "Sequence Evolution & Reconstruction"

// WARNING: this program is now deprecated. Please use 'reconstructSeq' as its replacement.
// Returns maximum likelihood sequence from an HBCGO primate alignment
package main

import (
	"flag"
	"fmt"
	"log"

	"github.com/vertgenlab/gonomics/dna"
	"github.com/vertgenlab/gonomics/exception"
	"github.com/vertgenlab/gonomics/expandedTree"
	"github.com/vertgenlab/gonomics/fasta"
	"github.com/vertgenlab/gonomics/reconstruct"
)

// likelihoodsToBaseUnbiased takes the un-normalized likelihoods for A, C, G, T as well
// and the probability
// threshold for when we will call it the mle base instead of N
// and gives back the reconstructed base for the hca.
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

// likelihoodsToBaseBias takes the un-normalized likelihoods for A, C, G, T as well
// as the index of the biased base (0 for A, 1 for C, etc), and the probability
// threshold for when we will call it the mle base instead of the biased base
// and gives back the reconstructed base for the hca/hga.
func likelihoodsToBaseBias(likes []float64, biasBase dna.Base, probThreshold float64, nonBiasProbThreshold float64) dna.Base {
	var total, nonHumanTotal, bestProb float64
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
			nonHumanTotal += likes[i]
		}
	}

	for i = range likes {
		if likes[i]/total >= probThreshold && likes[i] > bestProb && nonHumanTotal/total >= nonBiasProbThreshold {
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

func hcaIsPresent(human, bonobo, chimp, gorilla, orangutan dna.Base) bool {
	if baseIsPresent(human) && (baseIsPresent(bonobo) || baseIsPresent(chimp)) {
		return true
	}
	if (baseIsPresent(human) || baseIsPresent(bonobo) || baseIsPresent(chimp)) && (baseIsPresent(gorilla) || baseIsPresent(orangutan)) {
		return true
	}
	return false
}

func hgaIsPresent(human, bonobo, chimp, gorilla, orangutan dna.Base) bool {
	if baseIsPresent(gorilla) && (baseIsPresent(human) || baseIsPresent(chimp) || baseIsPresent(bonobo)) {
		return true
	}
	if baseIsPresent(orangutan) && (baseIsPresent(gorilla) || baseIsPresent(human) || baseIsPresent(chimp) || baseIsPresent(bonobo)) {
		return true
	}
	return false
}

func reconHcaBase(root, humanNode, chimpNode, nodeToRecon *expandedTree.ETree, position int, probThreshold float64, nonBiasProbThreshold float64, humanBias bool, chimpBias bool) {
	var likelihoods []float64
	var nextBase dna.Base
	reconstruct.SetState(root, position)
	likelihoods = reconstruct.FixFc(root, nodeToRecon)
	if humanBias {
		nextBase = likelihoodsToBaseBias(likelihoods, humanNode.Fasta.Seq[position], probThreshold, nonBiasProbThreshold)
	} else if chimpBias {
		nextBase = likelihoodsToBaseBias(likelihoods, chimpNode.Fasta.Seq[position], probThreshold, nonBiasProbThreshold)
	} else {
		nextBase = likelihoodsToBaseUnbiased(likelihoods, probThreshold)
	}
	nodeToRecon.Fasta.Seq = append(nodeToRecon.Fasta.Seq, nextBase)
}

func reconHgaBase(root, gorillaNode, nodeToRecon *expandedTree.ETree, position int, probThreshold float64, nonBiasProbThreshold float64) {
	var likelihoods []float64
	var nextBase dna.Base
	reconstruct.SetState(root, position)
	likelihoods = reconstruct.FixFc(root, nodeToRecon)
	nextBase = likelihoodsToBaseBias(likelihoods, gorillaNode.Fasta.Seq[position], probThreshold, nonBiasProbThreshold)
	nodeToRecon.Fasta.Seq = append(nodeToRecon.Fasta.Seq, nextBase)
}

func primateReconHcaMle(inFastaFilename string, inTreeFilename string, humanBias bool, chimpBias bool, probThreshold float64, nonHumanProbThreshold float64, useGenericNames bool, outputFastaFilename string) {
	var tree, humanNode, humanAltNode, bonoboNode, chimpNode, gorillaNode, orangutanNode, hcaNode *expandedTree.ETree
	var err error
	var i int

	tree, err = expandedTree.ReadTree(inTreeFilename, inFastaFilename)
	exception.FatalOnErr(err)
	if useGenericNames {
		humanNode = expandedTree.FindNodeName(tree, "human")
		if humanNode == nil {
			log.Fatalf("Didn't find human in the tree\n")
		}
		bonoboNode = expandedTree.FindNodeName(tree, "bonobo")
		if bonoboNode == nil {
			log.Fatalf("Didn't find bonobo in the tree\n")
		}
		chimpNode = expandedTree.FindNodeName(tree, "chimp")
		if chimpNode == nil {
			log.Fatalf("Didn't find chimp in the tree\n")
		}
		gorillaNode = expandedTree.FindNodeName(tree, "gorilla")
		if gorillaNode == nil {
			log.Fatalf("Didn't find gorilla in the tree\n")
		}
		orangutanNode = expandedTree.FindNodeName(tree, "orangutan")
		if orangutanNode == nil {
			log.Fatalf("Didn't find orangutan in the tree\n")
		}
	} else {
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
	}
	hcaNode = expandedTree.FindNodeName(tree, "hca")
	if hcaNode == nil {
		log.Fatalf("Didn't find hca in the tree\n")
	}

	if !(humanBias || chimpBias) {
		humanAltNode = expandedTree.FindNodeName(tree, "hg38alt")
		if humanAltNode == nil {
			log.Fatalf("Didn't find hg38alt in the tree\n")
		}
	}

	for i = range humanNode.Fasta.Seq {
		if hcaIsPresent(humanNode.Fasta.Seq[i], bonoboNode.Fasta.Seq[i], chimpNode.Fasta.Seq[i], gorillaNode.Fasta.Seq[i], orangutanNode.Fasta.Seq[i]) {
			reconHcaBase(tree, humanNode, chimpNode, hcaNode, i, probThreshold, nonHumanProbThreshold, humanBias, chimpBias)
		} else {
			hcaNode.Fasta.Seq = append(hcaNode.Fasta.Seq, dna.Gap)
		}
	}
	if humanBias || chimpBias {
		fasta.Write(outputFastaFilename, []fasta.Fasta{*humanNode.Fasta, *chimpNode.Fasta, *bonoboNode.Fasta, *gorillaNode.Fasta, *orangutanNode.Fasta, *hcaNode.Fasta})
	} else {
		fasta.Write(outputFastaFilename, []fasta.Fasta{*humanNode.Fasta, *humanAltNode.Fasta, *chimpNode.Fasta, *bonoboNode.Fasta, *gorillaNode.Fasta, *orangutanNode.Fasta, *hcaNode.Fasta})
	}
}

func primateReconHgaMle(inFile string, inTreeFileName string, probThreshold float64, nonBiasProbThreshold float64, useGenericNames bool, outFile string) {
	var tree, humanNode, bonoboNode, chimpNode, gorillaNode, orangutanNode, hgaNode *expandedTree.ETree
	var err error
	var i int

	tree, err = expandedTree.ReadTree(inTreeFileName, inFile)
	exception.FatalOnErr(err)
	if useGenericNames {
		humanNode = expandedTree.FindNodeName(tree, "human")
		if humanNode == nil {
			log.Fatalf("Didn't find human in the tree\n")
		}
		bonoboNode = expandedTree.FindNodeName(tree, "bonobo")
		if bonoboNode == nil {
			log.Fatalf("Didn't find bonobo in the tree\n")
		}
		chimpNode = expandedTree.FindNodeName(tree, "chimp")
		if chimpNode == nil {
			log.Fatalf("Didn't find chimp in the tree\n")
		}
		gorillaNode = expandedTree.FindNodeName(tree, "gorilla")
		if gorillaNode == nil {
			log.Fatalf("Didn't find gorilla in the tree\n")
		}
		orangutanNode = expandedTree.FindNodeName(tree, "orangutan")
		if orangutanNode == nil {
			log.Fatalf("Didn't find orangutan in the tree\n")
		}
	} else {
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
	}
	hgaNode = expandedTree.FindNodeName(tree, "hga")
	if hgaNode == nil {
		log.Fatalf("Didn't find hca in the tree\n")
	}

	for i = range humanNode.Fasta.Seq {
		if hgaIsPresent(humanNode.Fasta.Seq[i], bonoboNode.Fasta.Seq[i], chimpNode.Fasta.Seq[i], gorillaNode.Fasta.Seq[i], orangutanNode.Fasta.Seq[i]) {
			reconHgaBase(tree, gorillaNode, hgaNode, i, probThreshold, nonBiasProbThreshold)
		} else {
			hgaNode.Fasta.Seq = append(hgaNode.Fasta.Seq, dna.Gap)
		}
	}
	fasta.Write(outFile, []fasta.Fasta{*humanNode.Fasta, *chimpNode.Fasta, *bonoboNode.Fasta, *gorillaNode.Fasta, *orangutanNode.Fasta, *hgaNode.Fasta})
}

func primateRecon(infile string, outfile string, messyToN bool) {
	records := fasta.Read(infile)
	output := append(records, ParsimonyPrimateRecon(records, messyToN))
	fasta.Write(outfile, output)
}

// ParsimonyPrimateRecon returns a human-biased conservative human-chimp ancestor estimate from an input multiFa alignment in the order Human-Chimp-Bonoboo-Orangutan-Gorilla. MessyToN, when true, turns ambiguous bases to Ns in the output.
func ParsimonyPrimateRecon(records []fasta.Fasta, messyToN bool) fasta.Fasta {
	answer := fasta.Fasta{Name: "Human_Chimp_Ancestor"}

	var outputBase dna.Base

	if len(records) != 5 {
		log.Fatalf("Wrong number of sequences, expecting five, found %d.\n", len(records))
	}

	//confirm alignment lengths are all the same
	firstLength := len(records[0].Seq)
	for i := 1; i < len(records); i++ {
		if len(records[i].Seq) != firstLength {
			log.Fatalf("Sequence %d is the wrong length.\n", i+1)
		}
	}

	var Ncount int = 0
	var allMatch int = 0
	var humanMatchChimporBonobo int = 0
	var humanChange int = 0
	var gorillaVote int = 0
	var orangutanVote int = 0
	var messyBase int = 0

	human := records[0]
	bonobo := records[1]
	chimp := records[2]
	orangutan := records[3]
	gorilla := records[4]

	for j := 0; j < len(human.Seq); j++ {
		outputBase = human.Seq[j]
		if human.Seq[j] == dna.N { //output seq is N if human is N
			Ncount++
			outputBase = human.Seq[j]
		} else if humanInsertion(records, j) { //output seq is gap if human is insertion(see helper function)
			outputBase = dna.Gap
		} else if human.Seq[j] != dna.Gap && (chimp.Seq[j] == dna.Gap && bonobo.Seq[j] == dna.Gap) { //if there is sequence in humans and either gorilla or orangutan, but no sequence in chimp and bonobo, we output N, because we don't want to estimate the longer human to gorilla branch.
			if messyToN {
				outputBase = dna.N
			} else {
				outputBase = human.Seq[j]
			}
		} else if gorilla.Seq[j] == dna.Gap && orangutan.Seq[j] == dna.Gap { //no sequence outside the HCA (no gorilla or orangutan, ancestral state can't be determined)
			if messyToN {
				outputBase = dna.N
			} else {
				outputBase = human.Seq[j]
			}
		} else if human.Seq[j] == chimp.Seq[j] && human.Seq[j] == bonobo.Seq[j] { //if human matches chimp and bonobo, output is human
			outputBase = human.Seq[j]
			allMatch++
		} else if (human.Seq[j] == chimp.Seq[j] || human.Seq[j] == bonobo.Seq[j]) && human.Seq[j] != dna.Gap { //if human is a base and matches chimp or bonobo, output is human
			outputBase = human.Seq[j]
			humanMatchChimporBonobo++
		} else if (chimp.Seq[j] == bonobo.Seq[j] && (chimp.Seq[j] == gorilla.Seq[j] || chimp.Seq[j] == orangutan.Seq[j])) && (chimp.Seq[j] != dna.N && chimp.Seq[j] != dna.Gap) { //human is different from ancestor, chimp and bonobo agree with ancestor.
			outputBase = chimp.Seq[j]
			humanChange++
		} else if (human.Seq[j] == gorilla.Seq[j] || chimp.Seq[j] == gorilla.Seq[j] || bonobo.Seq[j] == gorilla.Seq[j]) && (gorilla.Seq[j] != dna.N && gorilla.Seq[j] != dna.Gap) { //more disagreement, but gorilla defines the ancestor in this case
			outputBase = gorilla.Seq[j]
			gorillaVote++
		} else if (human.Seq[j] == orangutan.Seq[j] || chimp.Seq[j] == orangutan.Seq[j] || bonobo.Seq[j] == orangutan.Seq[j] || gorilla.Seq[j] == orangutan.Seq[j]) && (orangutan.Seq[j] != dna.N && orangutan.Seq[j] != dna.Gap) {
			outputBase = orangutan.Seq[j]
			orangutanVote++
		} else {
			if human.Seq[j] != dna.Gap && !messyToN {
				outputBase = human.Seq[j]
			} else {
				outputBase = dna.N
			}
			messyBase++
		}
		answer.Seq = append(answer.Seq, outputBase)
	}
	return answer
}

func humanInsertion(records []fasta.Fasta, j int) bool {
	//true if the human sequence is the only one present at a position
	human := records[0]
	bonobo := records[1]
	chimp := records[2]
	orangutan := records[3]
	gorilla := records[4]

	if (human.Seq[j] != dna.Gap) && (chimp.Seq[j] == dna.Gap && bonobo.Seq[j] == dna.Gap && gorilla.Seq[j] == dna.Gap && orangutan.Seq[j] == dna.Gap) {
		return true
	}
	return false
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
	var mleHcaUnbiased *bool = flag.Bool("mleHcaUnbiased", false, "Estimates the human-chimpanzee ancestor. Default is an N, unless a base passes the probThreshold threshold.")
	var mleHcaHumanBiased *bool = flag.Bool("mleHcaHumanBiased", false, "Estimates the human-chimpanzee ancestor. Default is the human base, unless the non-human bases, collectively and individually, pass nonBiasProbThreshold and probThreshold.")
	var mleHcaChimpBiased *bool = flag.Bool("mleHcaChimpBiased", false, "Estimates the human-chimpanzee ancestor. Default is the chimp base, unless the non-chimp bases, collectively and individually, pass nonBiasPropThreshold and probThreshold.")
	var mleHgaGorillaBiased *bool = flag.Bool("mleHgaGorillaBiased", false, "Estimates the human-gorilla ancestor. Default is the gorilla base, unless the non-gorilla bases, collectively and individually, pass nonBiasPropThreshold and probThreshold.")
	var tree *string = flag.String("mle", "", "Filename for newick tree with branch lengths.  Must have the anticipated assembly names and the hca.")
	var probThreshold *float64 = flag.Float64("probThreshold", 0.0, "The probability that a base other than human must pass to be considered a true change in the hca.")
	var nonBiasProbThreshold *float64 = flag.Float64("nonBiasProbThreshold", 0.0, "The summation of all bases (other than the biased species) must pass this threshold for a non-biased species base to be considered in the hca.")
	var useGenericNames *bool = flag.Bool("useGenericNames", false, "Use generic names (human, bonobo, chimp, gorilla, orangutan) instead of assembly names (hg38, panPan2, panTro6, gorGor5, ponAbe3)")
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
	if *mleHcaHumanBiased && *mleHcaChimpBiased {
		log.Fatalf("Error: cannot be biased for both the human and the chimp base\n")
	}
	if *messyToN && (*mleHcaUnbiased || *mleHcaHumanBiased || *mleHcaChimpBiased) {
		log.Fatalf("Error: -messyToN can not be used with mle estimates\n")
	}
	if (*tree == "") && (*mleHcaUnbiased || *mleHcaHumanBiased || *mleHcaChimpBiased) {
		log.Fatalf("Error: you need to provide a tree when using an mle estimate\n")
	}
	if *mleHcaUnbiased && (*mleHcaHumanBiased || *mleHcaChimpBiased) {
		log.Fatalf("Error: Can not do both a biased and unbiased mle estimate\n")
	}
	if (*probThreshold != 0 || *nonBiasProbThreshold != 0) && !(*mleHcaUnbiased || *mleHcaHumanBiased || *mleHcaChimpBiased || *mleHgaGorillaBiased) {
		log.Fatalf("Error: Can not use probability threshold flags without also using an mle estimate\n")
	}
	if *nonBiasProbThreshold != 0 && *mleHcaUnbiased {
		log.Fatalf("Error: Can not do a nonBiasProbThreshold when also doing an unbiased estimate\n")
	}
	if *mleHgaGorillaBiased && (*mleHcaUnbiased || *mleHcaHumanBiased || *mleHcaChimpBiased) {
		log.Fatalf("Error: cannot estimate both the HCA and the HGA at once\n")
	}

	if *mleHcaUnbiased || *mleHcaHumanBiased || *mleHcaChimpBiased {
		// at this point we know that xor of the mleFlags is true, so we only pass mleHumanBiased and mleChimpBiased
		primateReconHcaMle(inFile, *tree, *mleHcaHumanBiased, *mleHcaChimpBiased, *probThreshold, *nonBiasProbThreshold, *useGenericNames, outFile)
	} else if *mleHgaGorillaBiased {
		primateReconHgaMle(inFile, *tree, *probThreshold, *nonBiasProbThreshold, *useGenericNames, outFile)
	} else {
		primateRecon(inFile, outFile, *messyToN)
	}
}
