package sam

import (
	"bytes"
	"encoding/hex"
	"io"
	"log"
	"math"
	"strconv"
	"strings"

	"github.com/vertgenlab/gonomics/bgzf"
	"github.com/vertgenlab/gonomics/cigar"
	"github.com/vertgenlab/gonomics/exception"
)

// BamWriter wraps a bgzf.BlockWriter and provides functions to write
// a Sam struct to a binary Bam file, including WriteToBamFileHandle.
type BamWriter struct {
	*bgzf.BlockWriter
	buf       bytes.Buffer // for storage between records
	recordBuf bytes.Buffer // for use within 1 record
	refMap    map[string]int
	u32       [4]byte // bytes for writing little endian
}

const (
	nul byte = '\u0000'
)

// NewBamWriter creates a new bam writer to the input io.Writer.
// The magic bam bytes and header are immediately written to the BamWriter.
func NewBamWriter(w io.Writer, h Header) *BamWriter {
	var bw BamWriter
	bw.BlockWriter = bgzf.NewBlockWriter(w)
	bw.buf.WriteString(magicBam)

	// write len of header text
	headerText := strings.Join(h.Text, "\n") + "\n"
	le.PutUint32(bw.u32[:], uint32(len(headerText)))
	bw.buf.Write(bw.u32[:])

	// write header text
	bw.buf.WriteString(headerText)

	// write number of ref sequences
	le.PutUint32(bw.u32[:], uint32(len(h.Chroms)))
	bw.buf.Write(bw.u32[:])

	// write reference information
	bw.refMap = make(map[string]int)
	for i, ref := range h.Chroms {
		// save order to refMap
		bw.refMap[ref.Name] = i

		// write len of ref name + 1 for nul byte
		le.PutUint32(bw.u32[:], uint32(len(ref.Name)+1))
		bw.buf.Write(bw.u32[:])

		// write ref name + nul byte
		bw.buf.WriteString(ref.Name)
		bw.buf.WriteByte(nul)

		// write ref len
		le.PutUint32(bw.u32[:], uint32(ref.Size))
		bw.buf.Write(bw.u32[:])
	}
	return &bw
}

// Close writes any data remaining in the buffer and closes
// the underlying bgzf.BlockWriter.
func (bw *BamWriter) Close() error {
	bytesRemain := bw.buf.Len()
	n, err := bw.Write(bw.buf.Bytes())
	if n != bytesRemain || err != nil {
		log.Panic(err)
	}
	return bw.BlockWriter.Close()
}

// WriteToBamFileHandle writes a single Sam struct to a bam file.
func WriteToBamFileHandle(bw *BamWriter, s Sam, bin uint16) {
	bw.recordBuf.Reset()

	// refID
	idx, ok := bw.refMap[s.RName]
	if s.RName == "*" {
		idx = -1
	} else if !ok {
		log.Fatalf("ERROR (WriteToBamFileHandle): Ref name '%s' "+
			"was present in read as RName, but not in header.", s.RName)
	}

	le.PutUint32(bw.u32[:4], uint32(idx))
	bw.recordBuf.Write(bw.u32[:4])

	// pos
	le.PutUint32(bw.u32[:4], s.Pos-1)
	bw.recordBuf.Write(bw.u32[:4])

	// len read name
	bw.recordBuf.WriteByte(uint8(len(s.QName) + 1))

	// mapq
	bw.recordBuf.WriteByte(s.MapQ)

	// BAI index bin
	le.PutUint16(bw.u32[:2], bin)
	bw.recordBuf.Write(bw.u32[:2])

	// num cigar op
	if cigar.IsUnmapped(s.Cigar) {
		le.PutUint16(bw.u32[:2], 0)
	} else {
		le.PutUint16(bw.u32[:2], uint16(len(s.Cigar)))
	}
	bw.recordBuf.Write(bw.u32[:2])

	// flag
	le.PutUint16(bw.u32[:2], s.Flag)
	bw.recordBuf.Write(bw.u32[:2])

	// len seq
	le.PutUint32(bw.u32[:4], uint32(len(s.Seq)))
	bw.recordBuf.Write(bw.u32[:4])

	// next ref id
	switch s.RNext {
	case "=": // same as RName
		// idx does not change

	case "*": // unmapped
		idx = -1

	default:
		idx, ok = bw.refMap[s.RNext]
		if !ok {
			log.Fatalf("ERROR (WriteToBamFileHandle): Ref name '%s' "+
				"was present in read as RNext, but not in header.", s.RName)
		}
	}

	le.PutUint32(bw.u32[:4], uint32(idx))
	bw.recordBuf.Write(bw.u32[:4])

	// next pos
	le.PutUint32(bw.u32[:4], s.PNext-1)
	bw.recordBuf.Write(bw.u32[:4])

	// tlen
	le.PutUint32(bw.u32[:4], uint32(s.TLen))
	bw.recordBuf.Write(bw.u32[:4])

	// read name nul terminated
	bw.recordBuf.WriteString(s.QName)
	bw.recordBuf.WriteByte(nul)

	// cigar is mapped (i.e. length is greater than 0).
	if !cigar.IsUnmapped(s.Cigar) {
		for i := range s.Cigar {
			le.PutUint32(bw.u32[:4], cigar.ToUint32(s.Cigar[i]))
			bw.recordBuf.Write(bw.u32[:4])
		}
	}
	// seq
	var seqInt uint8
	for i := range s.Seq {
		if i%2 == 0 {
			seqInt = baseEncoder[s.Seq[i]]
		} else {
			seqInt = seqInt << 4
			seqInt = seqInt | baseEncoder[s.Seq[i]]
			bw.recordBuf.WriteByte(seqInt)
		}
	}
	if len(s.Seq)%2 != 0 { // handle odd seq len
		seqInt = seqInt << 4
		bw.recordBuf.WriteByte(seqInt)
	}

	// qual
	if s.Qual == "*" {
		for range s.Seq {
			bw.recordBuf.WriteByte(0xff)
		}
	} else {
		for i := range s.Qual {
			bw.recordBuf.WriteByte(s.Qual[i] - 33)
		}
	}

	// aux data
	writeExtra(bw, s)

	// write block size then write recordBuf
	blockSize := uint32(bw.recordBuf.Len())
	le.PutUint32(bw.u32[:4], blockSize)
	bw.buf.Write(bw.u32[:4])           // write block size
	bw.buf.Write(bw.recordBuf.Bytes()) // write record to buf

	// Write to file once a full block has been read
	if bw.buf.Len() >= 64000 { // Max block size is 64KB
		n, err := bw.Write(bw.buf.Next(64000))
		if n != 64000 || err != nil {
			log.Panic(err)
		}
	}
}

// base        : =  A  C  M  G  R  S  V  T  W  Y  H  K  D  B  N
// bam specs   : 0  1  2  3  4  5  6  7  8  9  10 11 12 13 14 15
// gonomics val: 11 0  1  -1 2  -1 -1 -1 3  -1 -1 -1 -1 -1 -1 4
// lowercase converted to uppercase and all non-spec bases converted to 'N'.
var baseEncoder = []uint8{1, 2, 4, 8, 15, 1, 2, 4, 8, 15, 15, 15, 15, 15, 15, 15}

// writeExtra writes the Extra field of a sam to a BamWriter.
func writeExtra(bw *BamWriter, s Sam) {
	// write from unparsed Extra if present
	if len(s.unparsedExtra) != 0 {
		bw.recordBuf.Write(s.unparsedExtra)
		return
	}

	// return if no extra data
	if len(s.Extra) == 0 {
		return
	}

	// else parse tags as string
	tagSets := strings.Split(s.Extra, "\t")
	var triplet []string
	for i := range tagSets {
		triplet = retrieveTriplet(tagSets[i])
		writeTriplet(bw, triplet)
	}
}

// retrieveTriplet splits a tag string into type, number, and value.
func retrieveTriplet(tag string) []string {
	comp := strings.Split(tag, ":")
	if len(comp) < 3 || len(comp[0]) != 2 || len(comp[1]) != 1 { // checks to make sure tag is formatted properly
		log.Panicf("malformed auxiliary data '%s'", tag)
	}
	if comp[1] == "B" {
		comp[1] = comp[1] + comp[2][0:1] // the B case has another value as part of the array that needs to be stored in the second part of the triplet
		comp[2] = comp[2][2:]            // remove the character that goes to the second part of the triplet and the subsequent comma
	}
	comp[2] = strings.Join(comp[2:], ":") // in case value of tag contains literal ':'
	comp = comp[:3]
	return comp
}

// writeTriplet accepts a []string of len 3 that follows the format
// []string{ 2-byte tag, value type, values} and writes the parsed
// data to the input BamWriter per the Sam specifications.
func writeTriplet(bw *BamWriter, triplet []string) {
	if len(triplet) != 3 {
		log.Panicf("malformed auxiliary data '%s'", strings.Join(triplet, ":"))
	}

	// write tag bytes
	tag := triplet[0]
	if len(tag) != 2 {
		log.Panicf("auxiliary data tag must be exactly 2 characters offender:'%s'", tag)
	}
	bw.recordBuf.WriteString(tag)

	typ := triplet[1]
	values := strings.Split(triplet[2], ",")
	if len(values) == 1 && values[0] == "" {
		values = nil
	}
	realTyp := typ
	if typ[0] == 'B' { // 'B' is array of values
		bw.recordBuf.WriteByte(typ[0])
		// if typ is B, then the true data type is encoded as the second byte in typ
		realTyp = typ[1:2]
		bw.recordBuf.WriteByte(realTyp[0])
		le.PutUint32(bw.u32[:], uint32(len(values)))
		bw.recordBuf.Write(bw.u32[:])
	} else {
		bw.recordBuf.WriteByte(realTyp[0])
	}

	switch realTyp[0] {
	case 'A': // single character
		bw.recordBuf.WriteByte(values[0][0])

	case 'c', 'C': // int8, may be array
		var val int
		var err error
		for i := range values {
			val, err = strconv.Atoi(values[i])
			exception.PanicOnErr(err)
			bw.recordBuf.WriteByte(uint8(val))
		}

	case 's', 'S': // int16, may be array
		var val int
		var err error
		for i := range values {
			val, err = strconv.Atoi(values[i])
			exception.PanicOnErr(err)
			le.PutUint16(bw.u32[:2], uint16(val))
			bw.recordBuf.Write(bw.u32[:2])
		}

	case 'i', 'I': // int32, may be array
		var val int
		var err error
		for i := range values {
			val, err = strconv.Atoi(values[i])
			exception.PanicOnErr(err)
			le.PutUint32(bw.u32[:4], uint32(val))
			bw.recordBuf.Write(bw.u32[:4])
		}

	case 'f': // float32, may be array
		var val float64
		var err error
		for i := range values {
			val, err = strconv.ParseFloat(values[i], 32)
			exception.PanicOnErr(err)
			le.PutUint32(bw.u32[:4], math.Float32bits(float32(val)))
			bw.recordBuf.Write(bw.u32[:4])
		}

	case 'Z': // string, nul terminated
		bw.recordBuf.WriteString(strings.Join(values, ","))
		bw.recordBuf.WriteByte(nul)

	case 'H': // hex, nul terminated
		val, err := hex.DecodeString(strings.Join(values, ""))
		exception.PanicOnErr(err)
		bw.recordBuf.Write(val)
		bw.recordBuf.WriteByte(nul)

	default:
		log.Panicf("unrecognized auxiliary data type '%s'", typ)
	}
}
