package sam

import (
	"github.com/vertgenlab/gonomics/bgzf"
	"github.com/vertgenlab/gonomics/chromInfo"
	"io"
	"log"
	"strings"
)

// magicBam is a 4 byte sequence at the start of a bam file
const magicBam string = "BAM\u0001"

// BamReader wraps a bgzf.Reader with a fully allocated bgzf.Block.
type BamReader struct {
	bgzf.Reader
	blk *bgzf.Block
}

// OpenBam initiates a bgzf reader and parses the header info
// of the input bam file. OpenBam returns a fully initialized
// BamReader allocated with a bgzf Block. The second return is
// a Header struct parsed from plain header text which is stored
// in the bam file. The final return of []ChromInfo is redundant.
// Mainly with the Header return, but may be used if the plain
// text header is missing ref information in the bam file.
func OpenBam(filename string) (BamReader, Header, []chromInfo.ChromInfo) {
	var r BamReader
	r.Reader = bgzf.NewReader(filename)
	r.blk = bgzf.NewBlock()
	err := r.ReadBlock(r.blk)

	if err != nil && err != io.EOF { // EOF handled downstream
		log.Panic(err)
	}

	if r.blk.Len() == 0 {
		log.Fatalf("bam file empty: '%s'", filename)
	}

	var h Header
	var refs []chromInfo.ChromInfo
	refs, h.Text = parseBamHeader(r)
	return r, ParseHeaderText(h), refs
}

// parseBamHeader parses all header information in a bam file
// and returns the ref data as ChromInfo structs and header text.
func parseBamHeader(r BamReader) ([]chromInfo.ChromInfo, []string) {
	// check for magic bytes
	if string(r.blk.Next(4)) != magicBam {
		log.Fatal("missing magic bytes, bam file may be malformed")
	}

	// parse header text
	textLen := le.Uint32(r.blk.Next(4))
	text := string(r.blk.Next(int(textLen)))

	// parse references
	numRefs := int(le.Uint32(r.blk.Next(4)))
	refs := make([]chromInfo.ChromInfo, numRefs)

	var refNameLen int
	for i := 0; i < numRefs; i++ {
		refs[i].Order = i
		refNameLen = int(le.Uint32(r.blk.Next(4)))
		refs[i].Name = string(r.blk.Next(refNameLen))
		refs[i].Size = int(le.Uint32(r.blk.Next(4)))
	}

	lines := strings.Split(text, "\n")

	// in case plain text is newline terminated
	if lines[len(lines)-1] == "" {
		lines = lines[:len(lines)-1]
	}

	return refs, lines
}

// DecodeBam decodes a single Sam from a bgzf Block. The input
// block should be the output of a bgzf.ReadBlock either intact
// or with bytes discarded as determined by a virtual offset
// from a bai file. Regardless, the first byte in the input
// block should be the first byte in an alignment record.
func DecodeBam(r bgzf.Block) Sam {
	var s Sam

	return s
}
