package bed

import (
	"github.com/vertgenlab/gonomics/dna"
	"github.com/vertgenlab/gonomics/fasta"
	"log"
)

// ToLower converts bases in a fasta sequence to lowercase in specified bed regions.
func ToLower(fa fasta.Fasta, regions []Bed) {
	seq := fa.Seq
	for currRange := range regions {
		// assumes that the bed is formatted properly, ChromEnd > ChromStart
		// if asking for a region that exceeds the length of the base sequence, log fatal error
		if regions[currRange].ChromEnd > len(seq) {
			log.Fatalf("Error: ran out of sequence. Asking to manipulate a bed region where the exclusive/open-ended ChromEnd (%d) > length of sequence (%d).\n", regions[currRange].ChromEnd, len(seq))
		}
		dna.RangeToLower(seq, regions[currRange].ChromStart, regions[currRange].ChromEnd)
	}
}
