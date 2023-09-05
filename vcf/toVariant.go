package vcf

import (
	"log"

	"github.com/vertgenlab/gonomics/dna"
	"github.com/vertgenlab/gonomics/variant"
)

// SmallVariantChans wraps channels for all of the small variants that
// are valid variant.Mutators and variant.Effectors. All channels except
// for Records sends a parsed version of a given Vcf record that stores
// the corresponding variant type. Each send on one of the variant channels
// must be paired with a send on the records channel. i.e. for each parsed
// variant sent, the original vcf record is the following send. This is
// useful as the Vcf records often hold metadata that is not stored in the
// variant struct. Even if the original Vcf record is not used, the Records
// channel must be read else the sending goroutine will be blocked.
type SmallVariantChans struct {
	Substitutions chan variant.Substitution
	Insertions    chan variant.Insertion
	Deletions     chan variant.Deletion
	Delins        chan variant.Delins
	Records       chan Vcf
}

// NewSmallVariantChans creates a series of channels to send the variant
// representation of a given vcf record. The input buffSize determines
// the channel buffer size for all on the variant channels.
func NewSmallVariantChans(buffSize int) SmallVariantChans {
	var answer SmallVariantChans
	answer.Substitutions = make(chan variant.Substitution, buffSize)
	answer.Insertions = make(chan variant.Insertion, buffSize)
	answer.Deletions = make(chan variant.Deletion, buffSize)
	answer.Delins = make(chan variant.Delins, buffSize)
	answer.Records = make(chan Vcf, buffSize)
	return answer
}

// ChanToVariants splits an incoming channel of VCFs to a channel of variant types.
// Must be run as a goroutine else the thread will deadlock. Skips any vcf records
// where the ref/alt fields contains any of (:>[.). These records cannot be parsed
// into simple variants.
//
// See documentation for SmallVariantChans for more information on variant channels.
func ChanToVariants(c <-chan Vcf, sendChans SmallVariantChans) {
	var refSeq, altSeq []dna.Base
	for v := range c {
		if !canParseToBases(v.Ref) {
			continue // cannot parse vcf record so skip
		}
		refSeq = dna.StringToBases(v.Ref)
		for altIdx := range v.Alt {
			if !canParseToBases(v.Alt[altIdx]) {
				continue // cannot parse this allele so skip
			}
			altSeq = dna.StringToBases(v.Alt[altIdx])
			sendVariant(sendChans, v.Chr, v.Pos, refSeq, altSeq)
			sendChans.Records <- v
		}
	}
	closeSmallVariantChans(sendChans)
}

// GoChanToVariants wraps ChanToVariants and handles
// variant channel creation and goroutine spawning.
func GoChanToVariants(c <-chan Vcf) SmallVariantChans {
	answer := NewSmallVariantChans(cap(c))
	go ChanToVariants(c, answer)
	return answer
}

// sendVariant determines if a single alt allele derived from a vcf record is a
// Substitution, Insertion, Deletion, or Delins; parses the data to the appropriate
// variant type, and sends on the respective channel.
func sendVariant(sendChans SmallVariantChans, chr string, pos int, refSeq []dna.Base, altSeq []dna.Base) {
	var matchingOffset int
	refSeq, altSeq, matchingOffset = trimMatchingBases(refSeq, altSeq)
	pos += matchingOffset
	pos -= 1 // from 1-base to 0-base

	switch {
	case len(refSeq) == 1 && len(altSeq) == 1: // Substitution
		sendChans.Substitutions <- variant.Substitution{Chr: chr, Pos: pos, Ref: refSeq[0], Alt: altSeq[0]}

	case len(refSeq) < len(altSeq) && len(refSeq) == 0: // Insertion
		// e.g. ATG -> ATCG
		// VCF: pos = 2; ref = T; alt = TC
		// Desired: pos = 2; insSeq = C
		sendChans.Insertions <- variant.Insertion{Chr: chr, Pos: pos, Seq: altSeq}

	case len(refSeq) > len(altSeq) && len(altSeq) == 0: // Deletion
		// e.g. ATG -> AG
		// VCF: pos = 1; ref = AT; alt = A
		// Desired: start = 1; end = 2
		sendChans.Deletions <- variant.Deletion{Chr: chr, Start: pos, End: pos + len(refSeq)}

	case len(refSeq) > 0 && len(altSeq) > 0: // Delins
		// e.g. ATG -> AC
		// VCF: pos = 1; ref = ATG; alt = AC
		// Desired: start = 1; end = 3; insSeq = C
		sendChans.Delins <- variant.Delins{Chr: chr, Start: pos, End: pos + len(refSeq), InsSeq: altSeq}

	default:
		log.Panicf("could not parse identity of following variant\nChr: %s\nPos: %d\nRef: %s\nAlt: %s",
			chr, pos, dna.BasesToString(refSeq), dna.BasesToString(altSeq))
	}
}

// trimMatchingBases removed all left-aligned matching bases in ref and alt fields.
// returns the trimmed slices and the Number of bases trimmed.
func trimMatchingBases(a, b []dna.Base) ([]dna.Base, []dna.Base, int) {
	var offset int
	for len(a) > 0 && len(b) > 0 {
		if a[0] == b[0] {
			a = a[1:]
			b = b[1:]
			offset++
		} else {
			break
		}
	}
	if len(a) == 0 && len(b) == 0 {
		log.Panicf("error, all bases match trimmed in\n%s\nand\n%s", dna.BasesToString(a), dna.BasesToString(b))
	}
	return a, b, offset
}

// canParseToBases determines if a given string from a vcf file can be parsed to a []dna.Base.
func canParseToBases(s string) bool {
	for i := range s {
		if !validBase(s[i]) {
			return false
		}
	}
	return true
}

// validBase checks if a byte can be parsed to a dna.Base.
func validBase(b byte) bool {
	switch b {
	case 'A':
		return true
	case 'C':
		return true
	case 'G':
		return true
	case 'T':
		return true
	case 'a':
		return true
	case 'c':
		return true
	case 'g':
		return true
	case 't':
		return true
	default:
		return false
	}
}

// closeSmallVariantChans closes all channels wrapped by SmallVariantChans.
func closeSmallVariantChans(c SmallVariantChans) {
	close(c.Substitutions)
	close(c.Insertions)
	close(c.Deletions)
	close(c.Delins)
	close(c.Records)
}
