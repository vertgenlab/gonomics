package sam

import (
	"bytes"
	"errors"
	"github.com/vertgenlab/gonomics/bgzf"
	"github.com/vertgenlab/gonomics/chromInfo"
	"github.com/vertgenlab/gonomics/cigar"
	"github.com/vertgenlab/gonomics/dna"
	"io"
	"log"
	"strings"
)

// bam is a binary version of sam compressed as a bgzf file
// specs can be found in the sam documentation at:
// https://raw.githubusercontent.com/samtools/hts-specs/master/SAMv1.pdf

// magicBam is a 4 byte sequence at the start of a bam file
const magicBam string = "BAM\u0001"

// BamReader wraps a bgzf.BlockReader with a fully allocated bgzf.Block.
type BamReader struct {
	zr           *bgzf.BlockReader
	blk          *bgzf.Block
	intermediate bytes.Buffer
	refs         []chromInfo.ChromInfo
	eof          bool
}

// Close the BamReader and all underlying io.Readers.
func (r *BamReader) Close() error {
	return r.zr.Close()
}

// next returns the next n bytes from the bgzf block.
// next handles reading new block if more bytes are
// requested than are in the currently store block.
func (r *BamReader) next(n int) []byte {
	// we have enough bytes in intermediate
	if r.intermediate.Len() >= n {
		return r.intermediate.Next(n)
	}

	// we need to go to block for extra bytes
	if r.blk.Len() >= n-r.intermediate.Len() {
		if r.intermediate.Len() == 0 {
			return r.blk.Next(n)
		} else {
			blkLen := r.intermediate.Len()
			return append(r.intermediate.Next(n), r.blk.Next(n-blkLen)...)
		}
	}

	// not enough bytes in the currently stored block
	// read from multiple blocks if needed
	for r.intermediate.Len()+r.blk.Len() < n {
		_, err := r.intermediate.ReadFrom(r.blk)
		if err != nil {
			log.Panic(err)
		}

		err = r.zr.ReadBlock(r.blk)
		if err == io.EOF {
			r.eof = true
			break
		}
	}

	// return with fewest appends possible
	if r.intermediate.Len() != 0 {
		intLen := r.intermediate.Len()
		return append(r.intermediate.Next(n), r.blk.Next(n-intLen)...)
	}

	return r.blk.Next(n)
}

// OpenBam initiates a bgzf reader and parses the header info
// of the input bam file. OpenBam returns a fully initialized
// BamReader allocated with a bgzf Block. The second return is
// a Header struct parsed from plain header text which is stored
// in the bam file.
func OpenBam(filename string) (*BamReader, Header) {
	r := new(BamReader)
	r.zr = bgzf.NewBlockReader(filename)
	r.blk = bgzf.NewBlock()
	err := r.zr.ReadBlock(r.blk)
	if err != nil && err != io.EOF { // EOF handled downstream
		log.Panic(err)
	}

	if r.blk.Len() == 0 {
		log.Fatalf("bam file empty: '%s'", filename)
	}

	var h Header
	r.refs, h.Text = parseBamHeader(r)
	return r, ParseHeaderText(h)
}

// parseBamHeader parses all header information in a bam file
// and returns the ref data as ChromInfo structs and header text.
func parseBamHeader(r *BamReader) ([]chromInfo.ChromInfo, []string) {
	// check for magic bytes
	if string(r.next(4)) != magicBam {
		log.Fatal("missing magic bytes, bam file may be malformed")
	}

	// parse header text
	textLen := le.Uint32(r.next(4)) // le is an alias for binary.LittleEndian
	text := string(r.next(int(textLen)))

	// parse references
	numRefs := int(le.Uint32(r.next(4)))
	refs := make([]chromInfo.ChromInfo, numRefs)

	var refNameLen int
	for i := 0; i < numRefs; i++ {
		refs[i].Order = i
		refNameLen = int(le.Uint32(r.next(4)))
		// trim null from name
		refs[i].Name = trimNulOrPanic(string(r.next(refNameLen)))
		refs[i].Size = int(le.Uint32(r.next(4)))
	}

	lines := strings.Split(text, "\n")

	// in case plain text is newline terminated
	if lines[len(lines)-1] == "" {
		lines = lines[:len(lines)-1]
	}

	return refs, lines
}

// size in bytes of static portions of a bam alignment
// excludes block size uint32.
var staticBamAlnSize int = 32

// DecodeBam decodes a single Sam from a bgzf Block and stores the
// decoded sam in to input *Sam. Note that this overwrites all fields
// in the input sam, therefore if any information must be saved between
// decode calls, the values should be copied (e.g. copy slices instead of :=)
//
// The first return from DecodeBam is the binId that the decoded
// Sam was derived from. See sam.Bai for more information about
// bins and bam indexing. For most use cases, this value can be ignored.
//
// DecodeBam will return ErrNonStdBase if the input record contains
// bases other than A, C, G, T, or N. The sam return will still be
// fully parsed and can be used downstream, however all unsupported
// bases will be set to dna.Nil.
//
// At the end of the file, DecodeBam will return an empty sam and io.EOF.
func DecodeBam(r *BamReader, s *Sam) (binId uint32, err error) {
	blkSizeBytes := r.next(4) // bytes are read before using them to trigger EOF check
	if r.eof {
		return 0, io.EOF
	}
	blkSize := int(le.Uint32(blkSizeBytes)) // le is an alias for binary.LittleEndian
	refIdx := int32(le.Uint32(r.next(4)))
	if refIdx != -1 {
		s.RName = r.refs[refIdx].Name
	}
	s.Pos = le.Uint32(r.next(4)) + 1 // sam is 1 based
	lenReadName := int(r.next(1)[0])
	s.MapQ = r.next(1)[0]
	binId = uint32(le.Uint16(r.next(2))) // cast to uint32 since bai is a uint32
	numCigarOps := int(le.Uint16(r.next(2)))
	s.Flag = le.Uint16(r.next(2))
	lenSeq := int(le.Uint32(r.next(4)))
	refIdx = int32(le.Uint32(r.next(4)))
	s.RNext = "*"
	if refIdx != -1 {
		s.RNext = r.refs[refIdx].Name
	}
	if s.RNext == s.RName {
		s.RNext = "="
	}
	s.PNext = le.Uint32(r.next(4)) + 1 // sam is 1 based
	s.TLen = int32(le.Uint32(r.next(4)))

	// ******
	// Unsafe String Conversion
	/*
		qnameBytes := r.next(lenReadName)
		if qnameBytes[len(qnameBytes)-1] != 0 {
			log.Panicf("NUL byte not detected in QNAME\nBytes:%v", qnameBytes)
		}
		s.QName = unsafeByteToString(s.QName, qnameBytes[:len(qnameBytes)-1])
	*/
	// ******

	// ******
	// Safe String Conversion
	s.QName = trimNulOrPanic(string(r.next(lenReadName))) // safe version
	// ******

	if cap(s.Cigar) >= numCigarOps {
		s.Cigar = s.Cigar[:numCigarOps]
	} else {
		s.Cigar = make([]cigar.Cigar, numCigarOps)
	}
	var cigint uint32
	for i := 0; i < numCigarOps; i++ {
		cigint = le.Uint32(r.next(4))
		s.Cigar[i].Op = cigLookup[cigint&0xf]
		s.Cigar[i].RunLength = int(cigint >> 4)
	}

	// handle case where we are using a recycled sam struct.
	// in this case we don't want to waste memory and set the cigar to nil
	// for unaligned, so we use cig[0].Op = '*' which is how it is done when
	// reading from a sam file.
	if numCigarOps == 0 {
		if cap(s.Cigar) >= 1 {
			s.Cigar = s.Cigar[:1]
		} else {
			s.Cigar = make([]cigar.Cigar, 1)
		}
		s.Cigar[0].Op = '*'
	}

	if cap(s.Seq) >= lenSeq {
		s.Seq = s.Seq[:lenSeq]
	} else {
		s.Seq = make([]dna.Base, lenSeq)
	}
	err = readSeq(r, s.Seq)

	qual := r.next(lenSeq)
	for i := range qual {
		qual[i] += 33 // ascii offset for printable characters
	}

	// ******
	// Unsafe String Conversion
	//s.Qual = unsafeByteToString(s.Qual, qual) // unsafe version
	//if s.Qual[0] == 0xff {
	//	s.Qual = "*"
	//}
	// ******

	// ******
	// Safe String Conversion
	s.Qual = string(qual) // TODO this is 1 alloc per read, should change to []byte and remove unsafe ref above
	if len(qual) > 0 && qual[0]-33 == 0xff {
		s.Qual = "*"
	}
	// ******

	// The sam.Extra field is not parsed here as it would require parsing tags to their value, then
	// casting that value to a string. That would be excessively wasteful, so instead the bytes are
	// saved in the unparsedExtra field of sam, such that a user may call a function to parse these
	// bytes if and only if they need access to the tag fields. This should also make bam reading a
	// fair bit faster.
	ex := r.next(blkSize - (staticBamAlnSize +
		lenReadName + (4 * numCigarOps) + (((lenSeq) + 1) / 2) + lenSeq)) // to get remaining bytes in alignment
	if cap(s.unparsedExtra) < len(ex) {
		s.unparsedExtra = make([]byte, len(ex))
	}
	s.unparsedExtra = s.unparsedExtra[:len(ex)]
	copy(s.unparsedExtra, ex)
	s.parsedExtra = nil
	s.parsedExtraIdx = nil
	s.parsedExtraTags = nil
	s.parsedExtraTypes = nil
	return
}

// integer to cigar look quick lookup
var cigLookup = []rune{'M', 'I', 'D', 'N', 'S', 'H', 'P', '=', 'X', '*'}

// baseDecoder is for 4-bit to dna.Base decoding.
// gonomics does not support all 16 options in bam.
// options are: =ACMGRSVTWYHKDBN
var baseDecoder = [16]dna.Base{dna.Nil, dna.A, dna.C, dna.Nil, dna.G, dna.Nil, dna.Nil, dna.Nil, dna.T, dna.Nil, dna.Nil, dna.Nil, dna.Nil, dna.Nil, dna.Nil, dna.N}
var ErrNonStdBase error = errors.New("sequence contains bases other than A,C,G,T,N. Other bases are not supported in gonomics")

// readSeq reads bytes from r to fill the len of s with bases.
// Data in r are expected to be encoded with 4 bits per base.
func readSeq(r *BamReader, s []dna.Base) error {
	if len(s) == 0 {
		return nil
	}
	var err error
	var b byte
	var i int
	// limit is < len(s) - 2 so we can handle extra bits from odd lengths
	for i = 0; i < len(s)-2; i += 2 {
		b = r.next(1)[0]
		s[i] = baseDecoder[b>>4]
		s[i+1] = baseDecoder[b&0xf]
		if s[i] == dna.Nil || s[i+1] == dna.Nil {
			err = ErrNonStdBase
		}
	}

	b = r.next(1)[0]
	s[i] = baseDecoder[b>>4]
	if len(s)%2 == 0 { // seq is even, add both bases
		s[i+1] = baseDecoder[b&0xf]
	} // else extra bits are discarded
	return err
}

// trimNulOrPanic removes the NUL byte from NUL terminated strings.
// If the string does not end in NUL, the function with panic.
func trimNulOrPanic(s string) string {
	if s[len(s)-1] != 0 {
		log.Panicf("string '%s' is not NUL terminated. bam may be malformed", s)
	}
	return strings.TrimRight(s, "\u0000")
}

// Uncomment to use unsafe string conversion
// unsafeByteToString mutates the input string s to the string formed by
// calling `string(toCopy)`. This makes s unsafe for use between calls to
// DecodeBam unless s is copied to a new variable.
//func unsafeByteToString(s string, toCopy []byte) string {
//	header := (*reflect.StringHeader)(unsafe.Pointer(&s)) // get the header from s
//	bytesHeader := &reflect.SliceHeader{                  // convert to a byte slice header
//		Data: header.Data,
//		Len:  header.Len,
//		Cap:  header.Len,
//	}
//	if bytesHeader.Cap < len(toCopy) { // make a new slice if existing one is not big enough
//		return string(toCopy)
//	}
//
//	b := *(*[]byte)(unsafe.Pointer(bytesHeader)) // convert byte slice header to a literal []byte
//	header.Len = len(toCopy)                     // reduce len of slice header to actual len of toCopy
//	b = b[:len(toCopy)]                          // reduce len of literal byte slice conversion to len of toCopy
//	copy(b, toCopy)                              // mutate s
//	return s
//}
