// Command Group: "Data Conversion"

package main

import (
	"flag"
	"fmt"
	"github.com/vertgenlab/gonomics/cigar"
	"github.com/vertgenlab/gonomics/dna"
	"github.com/vertgenlab/gonomics/exception"
	"github.com/vertgenlab/gonomics/fasta"
	"github.com/vertgenlab/gonomics/fileio"
	"github.com/vertgenlab/gonomics/numbers"
	"github.com/vertgenlab/gonomics/sam"
	"github.com/vertgenlab/gonomics/vcf"
	"log"
	"strings"
	//"os"
	//"math/rand"
)

func samConsensus(samFileName string, refFile string, outFile string, vcfFile string) {
	ref := fasta.Read(refFile)
	fasta.AllToUpper(ref)

	votingMatrix := make(map[string][]voteBase, len(ref))
	var i, k int
	for i := range ref {
		votingMatrix[ref[i].Name] = make([]voteBase, len(ref[i].Seq))
		for k := range ref[i].Seq {
			votingMatrix[ref[i].Name][k] = voteBase{0, 0, 0, 0, 0}
		}
	}

	samFile := fileio.EasyOpen(samFileName)
	var err error
	var done bool = false
	var RefIndex, SeqIndex int
	var currentSeq []dna.Base
	var aln sam.Sam

	sam.ReadHeader(samFile)

	for aln, done = sam.ReadNext(samFile); done != true; aln, done = sam.ReadNext(samFile) {
		if aln.Cigar[0].Op != '*' {
			SeqIndex = 0
			RefIndex = int(aln.Pos - 1)
			for i = 0; i < len(aln.Cigar); i++ {
				currentSeq = aln.Seq
				if aln.Cigar[i].Op == 'D' {
					RefIndex = RefIndex + aln.Cigar[i].RunLength
				} else if cigar.CigarConsumesReference(*aln.Cigar[i]) {
					for k = 0; k < int(aln.Cigar[i].RunLength); k++ {
						switch currentSeq[SeqIndex] {
						case dna.A:
							votingMatrix[aln.RName][RefIndex].A++
						case dna.T:
							votingMatrix[aln.RName][RefIndex].T++
						case dna.G:
							votingMatrix[aln.RName][RefIndex].G++
						case dna.C:
							votingMatrix[aln.RName][RefIndex].C++
							//case dna.Gap:
							// votingMatrix[aln.RName][RefIndex].Gap++
						}
						SeqIndex++
						RefIndex++
					}
				} else if aln.Cigar[i].Op != 'H' {
					SeqIndex = SeqIndex + aln.Cigar[i].RunLength
				}
			}
		}
	}
	err = samFile.Close()
	exception.PanicOnErr(err)

	outVcfFile := fileio.EasyCreate(vcfFile)

	fmt.Fprintf(outVcfFile, "%s\n", strings.Join(vcf.NewHeader(samFileName).Text, "\n"))
	var current dna.Base
	var maxList []dna.Base
	for i := range ref {
		for k := range ref[i].Seq {
			current = voter(votingMatrix[ref[i].Name][k], maxList)
			if current == dna.N && ref[i].Seq[k] != dna.N {
				current = dna.ToLower(ref[i].Seq[k])
			} else if current != ref[i].Seq[k] {
				fmt.Fprintf(outVcfFile, "%s\t%d\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n", ref[i].Name, int64(k+1), ".", dna.BaseToString(ref[i].Seq[k]), dna.BaseToString(current), ".", ".", ".",".")
			}
			ref[i].Seq[k] = current
		}
	}
	err = outVcfFile.Close()
	exception.PanicOnErr(err)

	fasta.Write(outFile, ref)
}

func voter(v voteBase, maxList []dna.Base) dna.Base {
	if v.A+v.G+v.T+v.C+v.Gap == 0 {
		return dna.N
	}

	var max int32 = v.A
	var outBase dna.Base = dna.A

	if v.G > max {
		max = v.G
		outBase = dna.G
	}

	if v.C > max {
		max = v.C
		outBase = dna.C
	}

	if v.T > max {
		max = v.T
		outBase = dna.T
	}

	if v.Gap > max {
		max = v.Gap
		outBase = dna.Gap
	}

	//now we check for ties
	var maxCount int32 = 0

	if v.A == max {
		maxCount++
		maxList = append(maxList, dna.A)
	}
	if v.C == max {
		maxCount++
		maxList = append(maxList, dna.C)
	}
	if v.T == max {
		maxCount++
		maxList = append(maxList, dna.T)
	}
	if v.G == max {
		maxCount++
		maxList = append(maxList, dna.G)
	}
	if v.Gap == max {
		maxCount++
		maxList = append(maxList, dna.Gap)
	}

	if maxCount < 2 {
		return outBase
	}

	return maxList[numbers.RandIntInRange(0, len(maxList))]
}

type voteBase struct {
	A   int32
	C   int32
	G   int32
	T   int32
	Gap int32
}

func usage() {
	fmt.Print(
		"samConsensus - Generates a fasta file and accompanying vcf from a sam over a reference sequence.\n" +
			"Uncovered sequences are converted to lowercase reference sequences.\n" +
			"Usage:\n" +
			"samConsensus individual.sam ref.fa output.fa outputVariantList.vcf\n" +
			"options:\n")
	flag.PrintDefaults()
}

func main() {
	var expectedNumArgs int = 4
	flag.Usage = usage
	log.SetFlags(log.Ldate | log.Ltime | log.Lshortfile)
	flag.Parse()

	if len(flag.Args()) != expectedNumArgs {
		flag.Usage()
		log.Fatalf("Error: expecting %d arguments, but got %d\n",
			expectedNumArgs, len(flag.Args()))
	}

	inFile := flag.Arg(0)
	refFile := flag.Arg(1)
	outFile := flag.Arg(2)
	vcfFile := flag.Arg(3)

	samConsensus(inFile, refFile, outFile, vcfFile)
}
