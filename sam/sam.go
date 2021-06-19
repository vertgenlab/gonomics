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

// MatePair wraps to paired alignments into a single struct.
type MatePair struct {
	Fwd Sam
	Rev Sam
}

// Sam stores fields of a sam record in the corresponding data type
// denoted in the sam file specifications. The cigar and sequence
// fields are further parsed into a more complex data structure to
// facilitate ease of use.
type Sam struct {
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

	// TODO make the parse function
	// unparsedExtra is the Extra bytes from a bam file. If unparsedExtra != nil then
	// Extra is empty (by default). unparsedExtra can be parsed to values to access tag
	// field. If the struct was read directly from a sam file, unparsedExtra is never
	// used, therefore tag info must be parsed from the Extra field.
	unparsedExtra []byte
}

// Header encodes the header of a sam file as both the raw header (Text),
// and semi-parsed fields (Metadata and Chroms).
type Header struct {
	Text     []string
	Metadata Metadata
	Chroms   []chromInfo.ChromInfo // tags SQ - SN and SQ - LN
}

// String converts Sam to a string to satisfy the fmt.Stringer interface.
func (s Sam) String() string {
	return ToString(s)
}

// ToString converts an Sam struct to a tab delimited string per file specs.
func ToString(aln Sam) string {
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
