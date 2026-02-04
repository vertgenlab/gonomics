package bed

import (
	"github.com/vertgenlab/gonomics/dna"
	"github.com/vertgenlab/gonomics/fasta"
	"log"
	"strings"
)

// ToLower converts bases in a fasta sequence to lowercase in specified bed regions, where the bed chrom name matches the fasta record name
func ToLower(records []fasta.Fasta, regions []Bed, ignoreExtraRegions bool) {
	recordsMap := fasta.ToMap(records)
	var recordName string
	var currRecordSeq []dna.Base
	var found bool
	for currRegion := range regions { // loop through each bed region
		recordName = regions[currRegion].Chrom // bed chrom name matches the fasta record name
		// only convert if that fasta record name exists in the fasta recordsMap
		currRecordSeq, found = recordsMap[recordName]
		if found {
			// assumes that the bed is formatted properly, ChromEnd > ChromStart
			// if asking for a region that exceeds the length of the base sequence, log fatal error
			if regions[currRegion].ChromEnd > len(currRecordSeq) {
				log.Fatalf("Error: ran out of sequence. Asking to manipulate a bed region where ChromEnd (%d) > length of sequence (%d).\n", regions[currRegion].ChromEnd, len(currRecordSeq))
			}
			dna.RangeToLower(currRecordSeq, regions[currRegion].ChromStart, regions[currRegion].ChromEnd)
		} else {
			if !ignoreExtraRegions {
				log.Fatalf("Error: bed region chrom (%s) was not found as a fasta record name in the input fasta file\n", recordName)
			}
		}
	}
}

// SegregatingSites takes in a multiFa alignment and returns a new alignment containing only the columns with segregating sites, along with a bed file of the positions of segregating sites in the reference species.
// The inputs are the multiFa alignment (e.g. fasta records for human and hca), the chromosome the multiFa alignment is on (e.g. chr1 of human and hca), and a refStart offset integer value if necessary (e.g. the multiFa alignment is for a HAQER on chr1:100-200, so refStart is 100, and calculated SNP coordinates need to be offset by 100)
func SegregatingSites(aln []fasta.Fasta, chrom string, refStart int) ([]fasta.Fasta, []Bed) {
	// define variables
	var answerFa []fasta.Fasta = fasta.EmptyCopy(aln)
	var bedName string
	var currentBed Bed
	var answerBed []Bed

	speciesSeq := make([]string, len(aln))
	var chromStartAlnPos, chromStartRefPos int
	lastAlnPosConverted := 0
	lastRefPosConverted := 0

	// loop through multiFa
	for i := 0; i < len(aln[0].Seq); i++ {
		if fasta.IsSegregating(aln, i) {
			// report multiFa, collect base sequence in each species in preparation for reporting bed
			for k := 0; k < len(aln); k++ {
				answerFa[k].Seq = append(answerFa[k].Seq, aln[k].Seq[i])
				speciesSeq[k] = dna.BaseToString(aln[k].Seq[i])
				bedName = strings.Join(speciesSeq, "_")
			}
			// report bed entry for that 1 base position
			chromStartAlnPos = i
			chromStartRefPos = fasta.AlnPosToRefPosCounter(aln[0], chromStartAlnPos, lastRefPosConverted, lastAlnPosConverted)
			lastAlnPosConverted = chromStartAlnPos
			lastRefPosConverted = chromStartRefPos
			currentBed = Bed{Chrom: chrom, ChromStart: refStart + chromStartRefPos, ChromEnd: refStart + chromStartRefPos + 1, Name: bedName, Score: refStart + chromStartAlnPos, FieldsInitialized: 5} // Name field is referenceSpeciesBase_querySpeciesBase, Score field is AlnPos
			answerBed = append(answerBed, currentBed)
		}
	}
	return answerFa, answerBed
}
