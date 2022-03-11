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
				p.CountF[dna.Gap] = 0
				p.CountR[dna.Gap] = 0
				consensusBase = maxBase(p)
			}
		}
		if p.CountF[consensusBase] == 0 && p.CountR[consensusBase] == 0 { // skip if no data present
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
	var max int = p.CountF[dna.A] + p.CountR[dna.A]
	tiedBases := make([]dna.Base, 1) // zero value is dna.A

	max, tiedBases = getMax(p, max, dna.C, tiedBases)
	max, tiedBases = getMax(p, max, dna.G, tiedBases)
	max, tiedBases = getMax(p, max, dna.T, tiedBases)
	max, tiedBases = getMax(p, max, dna.Gap, tiedBases)

	if len(tiedBases) == 1 {
		return tiedBases[0]
	}

	return tiedBases[numbers.RandIntInRange(0, len(tiedBases))]
}

func getMax(p sam.Pile, currMax int, testBase dna.Base, tiedBases []dna.Base) (int, []dna.Base) {
	var count int = p.CountF[testBase] + p.CountR[testBase]
	if count > currMax {
		tiedBases = tiedBases[:1] // reset tied bases
		tiedBases[0] = testBase
		return count, tiedBases
	}

	if count == currMax {
		tiedBases = append(tiedBases, testBase)
	}

	return currMax, tiedBases
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
