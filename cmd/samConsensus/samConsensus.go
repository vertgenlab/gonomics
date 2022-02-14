// Command Group: "Data Conversion"

package main

import (
	"flag"
	"fmt"
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

func samConsensus(samFileName string, refFile string, outFile string, vcfFile string, skipGappedConsensus bool) {
	ref := fasta.Read(refFile)
	var err error
	for i := range ref {
		dna.AllToLower(ref[i].Seq)
	}
	refMap := fasta.ToMap(ref)

	reads, header := sam.GoReadToChan(samFileName)
	piles := sam.GoPileup(reads, header, false, nil, nil)
	outVcfFile := fileio.EasyCreate(vcfFile)

	_, err = fmt.Fprintf(outVcfFile, "%s\n", strings.Join(vcf.NewHeader(samFileName).Text, "\n"))
	exception.PanicOnErr(err)

	var consensusBase, refBase dna.Base
	var refName string
	for p := range piles {
		consensusBase = maxBase(p)
		if consensusBase == dna.Gap { // original code had commented out gap case
			if skipGappedConsensus { // ignore gapped position
				continue
			} else { // ignore gapped reads
				p.Count[dna.Gap] = 0
				consensusBase = maxBase(p)
			}
		}
		if p.Count[consensusBase] == 0 { // skip if no data present
			continue
		}
		refName = header.Chroms[p.RefIdx].Name
		refBase = dna.ToUpper(refMap[refName][p.Pos-1])
		refMap[header.Chroms[p.RefIdx].Name][p.Pos-1] = consensusBase
		if refBase != consensusBase {
			_, err = fmt.Fprintf(outVcfFile, "%s\t%d\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n", refName, int64(p.Pos), ".", dna.BaseToString(refBase), dna.BaseToString(consensusBase), ".", ".", ".", ".")
			exception.PanicOnErr(err)
		}
	}

	err = outVcfFile.Close()
	exception.PanicOnErr(err)

	fasta.Write(outFile, ref)
}

func maxBase(p sam.Pile) dna.Base {
	var max int = p.Count[0]
	var outBase dna.Base = dna.A

	if p.Count[dna.G] > max {
		max = p.Count[dna.G]
		outBase = dna.G
	}

	if p.Count[dna.C] > max {
		max = p.Count[dna.C]
		outBase = dna.C
	}

	if p.Count[dna.T] > max {
		max = p.Count[dna.T]
		outBase = dna.T
	}

	if p.Count[dna.Gap] > max {
		max = p.Count[dna.Gap]
		outBase = dna.Gap
	}

	//now we check for ties
	candidateBases := make([]int, 0, 5)
	if p.Count[dna.A] == max {
		candidateBases = append(candidateBases, int(dna.A))
	}
	if p.Count[dna.C] == max {
		candidateBases = append(candidateBases, int(dna.C))
	}
	if p.Count[dna.T] == max {
		candidateBases = append(candidateBases, int(dna.T))
	}
	if p.Count[dna.G] == max {
		candidateBases = append(candidateBases, int(dna.G))
	}
	if p.Count[dna.Gap] == max {
		candidateBases = append(candidateBases, int(dna.Gap))
	}

	if len(candidateBases) == 1 {
		return outBase
	}

	return dna.Base(candidateBases[numbers.RandIntInRange(0, len(candidateBases))])
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
	var skipGappedConsensus *bool = flag.Bool("skipGappedConsensus", false, "skipGappedConsensus skips regions where the consensus indicated a gapped position."+
		"This reports the position it would if no data were present. When false, this ignored gapped reads, but still determines the consensus sequence of any ungapped reads at that position.")
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

	samConsensus(inFile, refFile, outFile, vcfFile, *skipGappedConsensus)
}
