package main

import (
	"flag"
	"fmt"
	"log"
	"math/rand"

	"github.com/vertgenlab/gonomics/common"
	"github.com/vertgenlab/gonomics/exception"
	"github.com/vertgenlab/gonomics/fileio"
	"github.com/vertgenlab/gonomics/vcf"
)

func proximityBlockVcf(inFile string, outFile string, distance int, setSeed int64) {
	var err error
	rand.Seed(setSeed)
	records, header := vcf.Read(inFile)
	rand.Shuffle(len(records), func(i, j int) { records[i], records[j] = records[j], records[i] })

	//make an output slice and put the first entry in it.
	var retainedVcfs []vcf.Vcf = make([]vcf.Vcf, 1)
	retainedVcfs[0] = records[0]

	for i := 1; i < len(records); i++ {
		if vcfPassesDistanceThreshold(retainedVcfs, records[i], distance) {
			retainedVcfs = append(retainedVcfs, records[i])
		}
	}

	out := fileio.EasyCreate(outFile)
	vcf.NewWriteHeader(out, header)
	for _, j := range retainedVcfs {
		vcf.WriteVcf(out, j)
	}
	err = out.Close()
	exception.PanicOnErr(err)
}

func vcfPassesDistanceThreshold(retainedVcfs []vcf.Vcf, i vcf.Vcf, distanceThreshold int) bool {
	for _, j := range retainedVcfs {
		if i.Chr == j.Chr {
			if i.Pos < j.Pos {
				if j.Pos-i.Pos < distanceThreshold {
					return false
				}
			} else {
				if i.Pos-j.Pos < distanceThreshold {
					return false
				}
			}
		}
	}

	return true
}

func usage() {
	fmt.Print(
		"proximityBlockVcf - Pseudorandomly selects variants" +
			"from an input VCF file and retains variants in the output" +
			"that do not fall within a user-specified distance" +
			"of variants already chosen. Output is returned in a shuffled order.\n" +
			"Usage:\n" +
			"proximityBlockVcf input.vcf output.vcf distance" +
			"options:\n")
	flag.PrintDefaults()
}

func main() {
	var setSeed *int64 = flag.Int64("setSeed", -1, "Use a specific seed for the RNG.")

	var expectedNumArgs int = 3
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
	distance := common.StringToInt(flag.Arg(2))
	proximityBlockVcf(inFile, outFile, distance, *setSeed)
}
