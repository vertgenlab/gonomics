package sam

import (
	"bytes"
	"github.com/vertgenlab/gonomics/bgzf"
	"io"
	"log"
	"strings"
)

type BamWriter struct {
	bgzf.Writer
	buf       bytes.Buffer // for storage between records
	recordBuf bytes.Buffer // for use within 1 record
	refMap    map[string]int
}

const (
	nul byte = 0
)

// NewBamWriter creates a new bam writer to the input io.Writer.
// The magic bam bytes and header are immediately written to the BamWriter.
func NewBamWriter(w io.Writer, h Header) *BamWriter {
	var bw BamWriter
	bw.Writer = bgzf.NewWriter(w)
	bw.buf.WriteString(magicBam)

	// write len of header text
	headerText := strings.Join(h.Text, "\n")
	var u32 [4]byte
	le.PutUint32(u32[:], uint32(len(headerText)))
	bw.buf.Write(u32[:])

	// write header text
	bw.buf.WriteString(headerText)

	// write number of ref sequences
	le.PutUint32(u32[:], uint32(len(h.Chroms)))
	bw.buf.Write(u32[:])

	// write reference information
	bw.refMap = make(map[string]int)
	for i, ref := range h.Chroms {
		// save order to refMap
		bw.refMap[ref.Name] = i

		// write len of ref name + 1 for nul byte
		le.PutUint32(u32[:], uint32(len(ref.Name)+1))
		bw.buf.Write(u32[:])

		// write ref name + nul byte
		bw.buf.WriteString(ref.Name)
		bw.buf.WriteByte(nul)

		// write ref len
		le.PutUint32(u32[:], uint32(ref.Size))
		bw.buf.Write(u32[:])
	}
	return &bw
}

// Close writes any data remaining in the buffer and closes
// the underlying bgzf.Writer.
func (bw *BamWriter) Close() error {
	bytesRemain := bw.buf.Len()
	n, err := bw.Write(bw.buf.Bytes())
	if n != bytesRemain || err != nil {
		log.Panic(err)
	}
	return bw.Writer.Close()
}

// WriteToBamFileHandle writes a single Sam struct to a bam file.
func WriteToBamFileHandle(bw *BamWriter, s Sam, bin uint16) {
	bw.recordBuf.Reset()
	var u32 [4]byte

	// refID
	idx, ok := bw.refMap[s.RName]
	if !ok {
		idx = -1
	}
	le.PutUint32(u32[:4], uint32(idx))
	bw.recordBuf.Write(u32[:4])

	// pos
	le.PutUint32(u32[:4], s.Pos-1)
	bw.recordBuf.Write(u32[:4])

	// len read name
	le.PutUint32(u32[:4], uint32(len(s.QName)+1))
	bw.recordBuf.Write(u32[:4])

	// mapq
	bw.recordBuf.WriteByte(s.MapQ)

	// BAI index bin
	le.PutUint16(u32[:2], bin)
	bw.recordBuf.Write(u32[:2])

	// num cigar op
	le.PutUint16(u32[:2], uint16(len(s.Cigar)))
	bw.recordBuf.Write(u32[:2])

	// flag
	le.PutUint16(u32[:2], s.Flag)
	bw.recordBuf.Write(u32[:2])

	// len seq
	le.PutUint32(u32[:4], uint32(len(s.Seq)))
	bw.recordBuf.Write(u32[:4])

	// next ref id
	idx, ok = bw.refMap[s.RNext]
	if !ok {
		idx = -1
	}
	le.PutUint32(u32[:4], uint32(idx))
	bw.recordBuf.Write(u32[:4])

	// next pos
	le.PutUint32(u32[:4], s.PNext)
	bw.recordBuf.Write(u32[:4])

	// tlen
	le.PutUint32(u32[:4], uint32(s.TLen))
	bw.recordBuf.Write(u32[:4])

	// read name nul terminated
	bw.recordBuf.WriteString(s.QName)
	bw.recordBuf.WriteByte(nul)

	// TODO PICKUP AT CIGAR

	blockSize := uint32(bw.recordBuf.Len())
	le.PutUint32(u32[:4], blockSize)
	bw.buf.Write(u32[:4])              // write block size
	bw.buf.Write(bw.recordBuf.Bytes()) // write record to buf

	// Write to file once a full block has been read
	if bw.buf.Len() >= 64000 { // Max block size is 64KB
		n, err := bw.Write(bw.buf.Next(64000))
		if n != 64000 || err != nil {
			log.Panic(err)
		}
	}
}
