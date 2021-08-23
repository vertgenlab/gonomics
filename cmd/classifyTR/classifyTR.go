package main

import (
	"flag"
	"fmt"
	"github.com/vertgenlab/gonomics/fileio"
	"log"
	"os"
)

type eStrData struct { //maybe put this struct into a package for use when working with tandem repeats?
	Chrom				string
	EnsemblGene			string
	RepeatStart			int
	RepeatEnd			int
	MotifForward		string
	LocationToGene		string
	PutativeRole		string
	DistanceTrToGene	int
}

func classifyTR(eStrFile string, geneFile string, outFile string) {
	out := fileio.EasyCreate(outFile)
	eSTR type eStrData := aFunctionToParseTheInFileToTheEStrDataStruct
	// open the geneFile file
	geneInfo := make(map[string]int)
	geneInfo[geneFile$0] = int{geneFile$2} // I want to fill the map where column 1 value from the geneFile input file is
	// matched with the third column from the same file.
	// this geneFile has a header... does that matter? Does input file need to have it removed?
	//line code put the geneFile into the geneInfo map.


	for i := 1; i < len(eSTR); i ++ { //start at 1 to skip the header line.
		line of code to find the matching map to the eSTR.EnsembleGene value
		eStrMedian int := round((eSTR.RepeatStart + eSTR.RepeatEnd) / 2)

		TSS pos := look at column2 and find corresponding map/key combo

		if TSS pos > repeat median value{
		append to file or new column6: "upstream"
			if (TSS pos - eSTR median <= 2, 000){
				append or new column7: "promoter"
			} else
				append or new column7: "enhancer"
		}else{
			append to or new column6: "downstream"
			if (eSTR median - TSS pos <= 200) {
				append or new column7: "promoter"
			}else{
		append or new column7: "enhancer"
	}
	}
	basesApart := abs(TSS median - eSTR median)
		append or new column8: basesApart

	}

}

two input files.
	1) CSV
		column1: EnsID
		column2: chr
		column3: TSS pos
	- put this into a map where we have [ID] (key) --> coords (value)

2) TR file.
		column1: chr
		column2: ensembl
		column3: repeat start
		column4: repeat end
		column5: str motif forward
	- Ill loop over every line of this file to do the math.
		- find eSTR median value
		1. look at column2 and find corresponding map/key combo
		2. if TSS pos > repeat median value
				- append to file or new column6: "upstream"
					if (TSS pos - eSTR median <= 2,000)
							append or new column7: "promoter"
					else
							append or new column7: "enhancer"
			else
				- append to or new column6: "downstream"
					if (eSTR median - TSS pos <= 200)
						append or new column7: "promoter"
					else
						append or new column7: "enhancer"
		3. basesApart = abs(TSS median - eSTR median)
				append or new column8: basesApart


	- I can write output to a file?
			hold in a struct? where the struct is built to have 8 columns.
				- column6: str upstream or downstream of TSS?
				- column7: putative regulatory role? promoter or enhancer?
				- column8: bases from TSS?
			id really like to append to the input file, but not sure how?
	- this file will have a header unless I remove it. Remove it or figure out how to handle a header file?
in order to do this math, take median location of the str.
		= coordinate end + coordinate start / 2 ---- take absolute value of this output.

func usage() {
	fmt.Print(
		"classifyTR - Compares a genes TSS and the median coordinates of a tandem repeat to categorize a tandem repeat as a putative promoter or enhancer\n" +
			"Requires two input files containing tandem repeat data (chr, ensembl gene connected with, repeat start, repeat end,\n " +
			"and repeat motif forward) and ensembl gene data (ensembl id, chr, and TSS pos).\n " +
			"outoput file will be an appended version of the input tandemRepeatInfo.tab file"
			"tandemRepeatInfo.tab geneInfo.tab out.txt\n")
}


func main() {
	var expectedNumArgs int = 3
	flag.Usage = usage
	flag.Parse()

	if len(flag.Args()) != expectedNumArgs {
		flag.Usage()
		log.Fatalf("Error: expecting %d arguments, but got %d\n",
			expectedNumArgs, len(flag.Args()))
	}

	eStrFile := flag.Args(0)
	geneFile := flag.Args(1)
	outFile := flag.Args(2)
	classifyTR(eStrFile, geneFile, outFile)


}
