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

// SegregatingSites takes in a multiFa alignment and returns a new alignment containing only the columns with segregating sites, along with a bed file of the positions of segregating sites in the reference species
func SegregatingSites(aln []fasta.Fasta) ([]fasta.Fasta, []int, []string) {
	// define variables
	var answer []fasta.Fasta = fasta.EmptyCopy(aln)
	var i, k int
	var bedPos []int
	speciesSeq := make([]string, len(aln))
	var bedName string
	var bedNames []string

	// loop through multiFa
	for i = 0; i < len(aln[0].Seq); i++ {
		if fasta.IsSegregating(aln, i) {
			// report multiFa, collect base sequence in each species in preparation for reporting bed
			for k = 0; k < len(aln); k++ {
				answer[k].Seq = append(answer[k].Seq, aln[k].Seq[i])
				speciesSeq[k] = dna.BaseToString(aln[k].Seq[i])
				bedName = strings.Join(speciesSeq, "_")
			}
			// report bed entry for that 1 base position
			bedPos = append(bedPos, i)
			bedNames = append(bedNames, bedName)
		}
	}
	return answer, bedPos, bedNames
}
