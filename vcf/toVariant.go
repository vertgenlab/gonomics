package vcf

import (
	"github.com/vertgenlab/gonomics/dna"
	"github.com/vertgenlab/gonomics/variant"
	"log"
	"strings"
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
// Must be run as a goroutine else the thread will deadlock.
// See documentation for SmallVariantChans for more information on variant channels.
func ChanToVariants(c <-chan Vcf, sendChans SmallVariantChans) {
	for v := range c {
		if !canParseToBases(v.Ref) {
			continue // cannot parse vcf record so skip
		}
		var refSeq, altSeq []dna.Base
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

func trimMatchingBases(a, b []dna.Base) (aNew []dna.Base, bNew []dna.Base, offset int) {
	aNew, bNew = a, b
	for len(aNew) > 0 && len(bNew) > 0 {
		if aNew[0] == bNew[0] {
			aNew = aNew[1:]
			bNew = bNew[1:]
			offset++
		} else {
			break
		}
	}
	if len(aNew) == 0 && len(bNew) == 0 {
		log.Panicf("error, all bases match trimmed in\n%s\nand\n%s", dna.BasesToString(a), dna.BasesToString(b))
	}
	return
}

// canParseToBases determines if a given string from a vcf file can be parsed to a []dna.Base
func canParseToBases(s string) bool {
	return !strings.ContainsAny(s, ":>[.")
}

// closeSmallVariantChans closes all channels wrapped by SmallVariantChans
func closeSmallVariantChans(c SmallVariantChans) {
	close(c.Substitutions)
	close(c.Insertions)
	close(c.Deletions)
	close(c.Delins)
	close(c.Records)
}
