package bed

import (
	"github.com/vertgenlab/gonomics/dna"
	"github.com/vertgenlab/gonomics/fasta"
	"log"
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
