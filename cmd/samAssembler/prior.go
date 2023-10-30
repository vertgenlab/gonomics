package main

import (
	"flag"
	"fmt"
	"github.com/vertgenlab/gonomics/dna"
	"github.com/vertgenlab/gonomics/exception"
	"github.com/vertgenlab/gonomics/fasta"
	"github.com/vertgenlab/gonomics/fileio"
	"github.com/vertgenlab/gonomics/sam"
)

type PriorSettings struct {
	SamFileName         string
	ReferenceFile       string
	OutFile             string
	Epsilon             float64
	LikelihoodCacheSize int
	PseudoCount         float64
	AsCounts            bool
}

func priorUsage() {
	fmt.Print("samAssembler prior - Construct an empirical conditional Dirichlet prior for output diploid" +
		"genotypes based on maximum likelihood estimation of genotypes from aligned short reads.\n" +
		"Usage: \n" +
		"samAssembler [options] prior reads.sam/bam ref.fa output.txt\n" +
		"Options:\n")
	flag.PrintDefaults()
}

func SamAssemblerPrior(s PriorSettings) {
	var answer = make([][]float64, 4)
	for i := range answer {
		answer[i] = make([]float64, 10)
		for j := range answer[i] {
			answer[i][j] = s.PseudoCount
		}
	}
	var refBase dna.Base
	var currChrom string
	var baseCall sam.DiploidBase
	var err error

	// read pileups from sam/bam
	reads, header := sam.GoReadToChan(s.SamFileName)
	piles := sam.GoPileup(reads, header, false, nil, nil)
	//read reference from fasta
	ref := fasta.Read(s.ReferenceFile)
	for i := range ref {
		dna.AllToUpper(ref[i].Seq)
	}
	refMap := fasta.ToMap(ref)

	homozygousCache := make([][]float64, s.LikelihoodCacheSize)
	for i := range homozygousCache {
		homozygousCache[i] = make([]float64, s.LikelihoodCacheSize)
	}
	heterozygousCache := make([][]float64, s.LikelihoodCacheSize)
	for i := range heterozygousCache {
		heterozygousCache[i] = make([]float64, s.LikelihoodCacheSize)
	}
	diploidBasePriorCache := sam.MakeDiploidBaseFlatPriorCache()

	for p := range piles {
		currChrom = header.Chroms[p.RefIdx].Name
		refBase = refMap[currChrom][p.Pos-1]
		if refBase < 4 {
			baseCall = sam.DiploidBaseCallFromPile(p, refBase, diploidBasePriorCache, homozygousCache, heterozygousCache, sam.AncientLikelihoodCache{}, s.Epsilon, 0)
			if baseCall < 10 {
				answer[refBase][baseCall]++
			}
		}
	}

	if !s.AsCounts {
		answer = convertToProb(answer)
	}

	out := fileio.EasyCreate(s.OutFile)
	_, err = fmt.Fprintf(out, ".\tAA\tAC\tAG\tAT\tCC\tCG\tCT\tGG\tGT\tTT\n")
	exception.PanicOnErr(err)
	_, err = fmt.Fprintf(out, "RefA\t%v\t%v\t%v\t%v\t%v\t%v\t%v\t%v\t%v\t%v\n",
		answer[dna.A][0],
		answer[dna.A][1],
		answer[dna.A][2],
		answer[dna.A][3],
		answer[dna.A][4],
		answer[dna.A][5],
		answer[dna.A][6],
		answer[dna.A][7],
		answer[dna.A][8],
		answer[dna.A][9],
	)
	exception.PanicOnErr(err)
	_, err = fmt.Fprintf(out, "RefC\t%v\t%v\t%v\t%v\t%v\t%v\t%v\t%v\t%v\t%v\n",
		answer[dna.C][0],
		answer[dna.C][1],
		answer[dna.C][2],
		answer[dna.C][3],
		answer[dna.C][4],
		answer[dna.C][5],
		answer[dna.C][6],
		answer[dna.C][7],
		answer[dna.C][8],
		answer[dna.C][9],
	)
	exception.PanicOnErr(err)
	_, err = fmt.Fprintf(out, "RefG\t%v\t%v\t%v\t%v\t%v\t%v\t%v\t%v\t%v\t%v\n",
		answer[dna.G][0],
		answer[dna.G][1],
		answer[dna.G][2],
		answer[dna.G][3],
		answer[dna.G][4],
		answer[dna.G][5],
		answer[dna.G][6],
		answer[dna.G][7],
		answer[dna.G][8],
		answer[dna.G][9],
	)
	exception.PanicOnErr(err)
	_, err = fmt.Fprintf(out, "RefT\t%v\t%v\t%v\t%v\t%v\t%v\t%v\t%v\t%v\t%v\n",
		answer[dna.T][0],
		answer[dna.T][1],
		answer[dna.T][2],
		answer[dna.T][3],
		answer[dna.T][4],
		answer[dna.T][5],
		answer[dna.T][6],
		answer[dna.T][7],
		answer[dna.T][8],
		answer[dna.T][9],
	)
	exception.PanicOnErr(err)

	err = out.Close()
	exception.PanicOnErr(err)
}

func convertToProb(input [][]float64) [][]float64 {
	var currRowSum float64
	var currRow, currColumn int
	var output [][]float64 = make([][]float64, len(input))
	for currRow = range output {
		output[currRow] = make([]float64, len(input[currRow]))
	}
	for currRow = range input {
		currRowSum = 0
		for currColumn = range input[currRow] {
			currRowSum += input[currRow][currColumn]
		}
		for currColumn = range input[currRow] {
			output[currRow][currColumn] = input[currRow][currColumn] / currRowSum
		}
	}
	return output
}
