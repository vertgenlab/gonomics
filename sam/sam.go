// Package sam provides functions for reading, writing, and manipulating sam files.
// The Sam struct is built off the v1.6 specification defined at https://github.com/samtools/hts-specs
package sam

import (
	"strconv"
	"strings"

	"github.com/vertgenlab/gonomics/chromInfo"
	"github.com/vertgenlab/gonomics/cigar"
	"github.com/vertgenlab/gonomics/dna"
	"github.com/vertgenlab/gonomics/exception"
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
	QName string        // query template name: [!-?A-~]{1,254}
	Flag  uint16        // bitwise flag, bits defined in info.go
	MapQ  uint8         // mapping quality
	RName string        // reference sequence name: \*|[:rname:∧ *=][:rname:]*
	Pos   uint32        // 1-based leftmost mapping position
	Cigar []cigar.Cigar // parsed cigar: originally string with \*|([0-9]+[MIDNSHPX=])+
	RNext string        // reference name of the mate/next read: \*|=|[:rname:∧ *=][:rname:]*
	PNext uint32        // position of the mate/next read
	TLen  int32         // observed template length
	Seq   []dna.Base    // parsed sequence: originally string with \*|[A-Za-z=.]+
	Qual  string        // ASCII of Phred-scaled base quality+33: [!-~]+ // TODO: parse to []Qual???
	Extra string        // TODO: parse to map or slice w/ index embedded in header???

	// unparsedExtra is the Extra bytes from a bam file. If unparsedExtra != nil then
	// Extra is empty (by default). unparsedExtra can be parsed to values to access tag
	// field. If the struct was read directly from a sam file, unparsedExtra is never
	// used, therefore tag info must be parsed from the Extra field.
	unparsedExtra []byte

	// parsedExtra stores the tag information keying each tag (e.g. RG, RX, NM) to that
	// tags value. The map returns an integer that corresponds to the index in parsedExtra,
	// which stores the underlying value of either: uint8, int8, uint16, int16, uint32,
	// int32, float32, or a string. parsedExtra is nil until the first QueryTag call which
	// will trigger the parsing of unparsedExtra.
	parsedExtraIdx   map[Tag]int
	parsedExtraTags  []Tag
	parsedExtraTypes [][]byte // first idx is idx in parsed extra, second index is for B is any of "cCsSiIfZH"
	parsedExtra      []interface{}
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
	var buf strings.Builder

	var err error
	if len(s.unparsedExtra) > 0 {
		s, err = parseExtra(s)
		exception.PanicOnErr(err)
		s.Extra = parsedExtraToString(&s)
	}
	buf.WriteString(s.QName)
	buf.WriteByte('\t')

	buf.WriteString(strconv.Itoa(int(s.Flag)))
	buf.WriteByte('\t')

	buf.WriteString(s.RName)
	buf.WriteByte('\t')

	buf.WriteString(strconv.Itoa(int(s.Pos)))
	buf.WriteByte('\t')

	buf.WriteString(strconv.Itoa(int(s.MapQ)))
	buf.WriteByte('\t')

	buf.WriteString(cigar.ToString(s.Cigar))
	buf.WriteByte('\t')

	buf.WriteString(s.RNext)
	buf.WriteByte('\t')
	buf.WriteString(strconv.Itoa(int(s.PNext)))
	buf.WriteByte('\t')

	buf.WriteString(strconv.Itoa(int(s.TLen)))
	buf.WriteByte('\t')

	buf.WriteString(dna.BasesToString(s.Seq))
	buf.WriteByte('\t')

	buf.WriteString(s.Qual)

	if s.Extra != "" {
		buf.WriteByte('\t')
		buf.WriteString(s.Extra)
	}

	return buf.String()
}

// 	buf.WriteString(s.Qual)

// 	if len(s.unparsedExtra) > 0 {
//         var err error
//         s, err = parseExtra(s)
//         if err != nil {
//             // Handle the error gracefully, perhaps by logging or using a default value
//             // For now, let's just skip the extra field
//             fmt.Printf("Error parsing extra fields: %v\n", err)
//         } else {
//             extraStr := parsedExtraToString(&s)
//             if extraStr != "" {
//                 buf.WriteByte('\t')
//                 buf.WriteString(extraStr)
//             }
//         }
//     }

// 	return buf.String()
// }
