package fastq

import (
	"fmt"

	"github.com/vertgenlab/gonomics/dna"
	"github.com/vertgenlab/gonomics/exception"
	"github.com/vertgenlab/gonomics/fileio"
)

// SingleCellPair is a struct that contains single-cell sequencing information, including a pairedEnd read, barcode, and UMI.
type SingleCellPair struct {
	Reads PairedEnd
	Bx    []dna.Base
	Umi   []dna.Base
}

// ReadToChanSingleCellPair takes two files, corresponding to R1 and R2, and returns a channel of SingleCellPair structs. It parses the barcode and UMI to user-specified lengths.
func ReadToChanSingleCellPair(fileOne string, fileTwo string, barcodeLength int, umiLength int, tenXg chan<- SingleCellPair) {
	fwd, rev := fileio.EasyOpen(fileOne), fileio.EasyOpen(fileTwo)

	for curr, done := NextFastqPair(fwd, rev); !done; curr, done = NextFastqPair(fwd, rev) {
		tenXg <- PairedEndToSingleCellPair(curr, barcodeLength, umiLength)
	}
	var err error
	err = fwd.Close()
	exception.PanicOnErr(err)
	err = rev.Close()
	exception.PanicOnErr(err)
	close(tenXg)
}

// PairedEndToSingleCellPair takes a PairedEnd read struct and parses a SingleCellPair with user-specified barcode and Umi lengths.
func PairedEndToSingleCellPair(fqPair PairedEnd, barcodeLength int, umiLength int) SingleCellPair {
	sc := SingleCellPair{Reads: fqPair, Bx: GetBarcode(fqPair.Fwd, 0, barcodeLength), Umi: GetBarcode(fqPair.Fwd, barcodeLength, barcodeLength+umiLength)}
	sc.Reads.Fwd = TrimFastq(fqPair.Fwd, barcodeLength+umiLength, len(fqPair.Fwd.Seq))
	bxTag := fmt.Sprintf("BX:%s", dna.BasesToString(sc.Bx))
	umiTag := fmt.Sprintf("UMI:%s", dna.BasesToString(sc.Umi))
	sc.Reads.Fwd.Name = fmt.Sprintf("%s_%s_%s", fqPair.Fwd.Name, umiTag, bxTag)
	sc.Reads.Rev.Name = fmt.Sprintf("%s_%s_%s", fqPair.Rev.Name, umiTag, bxTag)
	return sc
}

// TrimFastq trims fastq bases and quals given a start and end position zero base. Left-closed right-open.
func TrimFastq(fq Fastq, start int, end int) Fastq {
	fq.Seq = fq.Seq[start:end]
	fq.Qual = fq.Qual[start:end]
	return fq
}

// GetBarcode copies a subset of the given fastq record, at start and and end sites, zero base
// could be a UMI or 10x barcode.
func GetBarcode(fq Fastq, start int, end int) []dna.Base {
	answer := make([]dna.Base, end-start)
	copy(answer, fq.Seq[start:end])
	return answer
}

// UmiCollision takes a SingleCellPair and consults a map (with the string as the concatenated barcode and UMI strings), \\
// Returns true if the SingleCellPair UMI and barcode are found in the map, false otherwise. If a pair is not found in the map, the map is updated to include the current pair.
func UmiCollision(sc SingleCellPair, umiMap map[string]int) bool {
	currUmi := fmt.Sprintf("%s%s", sc.Bx, sc.Umi)
	count, found := umiMap[currUmi]
	if !found {
		umiMap[currUmi] = 1
	} else {
		umiMap[currUmi] = count + 1
	}
	return found
}
