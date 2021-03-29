// Package sam provides functions for reading, writing, and manipulating sam files.
// The Sam struct is built off the v1.6 specification defined at https://github.com/samtools/hts-specs
package sam

import (
	"fmt"
	"github.com/vertgenlab/gonomics/chromInfo"
	"github.com/vertgenlab/gonomics/cigar"
	"github.com/vertgenlab/gonomics/dna"
)

const (
	samSpecVersion = "1.6" // Version of sam this package is built on
)

type Aln struct {
	QName string         // query template name: [!-?A-~]{1,254}
	Flag  uint16         // bitwise flag, bits defined in info.go
	MapQ  uint8          // mapping quality
	RName string         // reference sequence name: \*|[:rname:∧ *=][:rname:]*
	Pos   uint32         // 1-based leftmost mapping position
	Cigar []*cigar.Cigar // parsed cigar: originally string with \*|([0-9]+[MIDNSHPX=])+
	RNext string         // reference name of the mate/next read: \*|=|[:rname:∧ *=][:rname:]*
	PNext uint32         // position of the mate/next read
	TLen  int32          // observed template length
	Seq   []dna.Base     // parsed sequence: originally string with \*|[A-Za-z=.]+
	Qual  string         // ASCII of Phred-scaled base quality+33: [!-~]+ // TODO: parse to []Qual???
	Extra string         // TODO: parse to map or slice w/ index embedded in header???
}

// Header encodes the header of a sam file as both the raw header (Text),
// and semi-parsed fields (Metadata and Chroms).
type Header struct {
	Text     []string
	Metadata Metadata
	Chroms   []chromInfo.ChromInfo // tags SQ - SN and SQ - LN
}

// ToString converts an Aln struct to a tab delimited string per file specs.
func ToString(aln Aln) string {
	var answer string
	if aln.Extra == "" {
		answer = fmt.Sprintf("%s\t%d\t%s\t%d\t%d\t%s\t%s\t%d\t%d\t%s\t%s",
			aln.QName, aln.Flag, aln.RName, aln.Pos, aln.MapQ, cigar.ToString(aln.Cigar), aln.RNext, aln.PNext, aln.TLen, dna.BasesToString(aln.Seq), aln.Qual)
	} else {
		answer = fmt.Sprintf("%s\t%d\t%s\t%d\t%d\t%s\t%s\t%d\t%d\t%s\t%s\t%s",
			aln.QName, aln.Flag, aln.RName, aln.Pos, aln.MapQ, cigar.ToString(aln.Cigar), aln.RNext, aln.PNext, aln.TLen, dna.BasesToString(aln.Seq), aln.Qual, aln.Extra)
	}
	return answer
}
