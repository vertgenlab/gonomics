package sam

import (
	"github.com/vertgenlab/gonomics/chromInfo"
	"github.com/vertgenlab/gonomics/exception"
	"io"
	"log"
)

// SeekBamRegion returns a slice of reads that overlap the input region. SeekBamRegion will advance the
// reader as necessary so beware continuing linear reading after a call to SeekBamRegion.
func SeekBamRegion(br *BamReader, bai Bai, chrom string, start, end uint32) []Sam {
	if start > end {
		log.Panicf("ERROR: SeekBamRegion input start > end. %d > %d\n", start, end)
	}
	var err error
	var coffset, uoffset uint64
	var ans []Sam
	refIdx := chromInfo.SliceToMap(br.refs)[chrom].Order
	bins := regionToBins(int(start), int(end))
	ref := bai.refs[refIdx]
	var i, j int
	var ok bool
	for i = range bins { // retrieve bins that may contain overlapping reads
		if _, ok = ref.binIdIdx[bins[i]]; !ok {
			continue
		}
		for j = range ref.bins[ref.binIdIdx[i]].chunks { // check all chunks with each bin
			uoffset = ref.bins[ref.binIdIdx[i]].chunks[j].start & 0xFFFF // byte offset into uncompressed data stream
			coffset = ref.bins[ref.binIdIdx[i]].chunks[j].start >> 16    // byte offset from start to bgzf block
			_, err = br.zr.Seek(int64(coffset), io.SeekStart)
			exception.PanicOnErr(err)
			err = br.zr.ReadBlock(br.blk)
			exception.PanicOnErr(err)
			br.blk.Next(int(uoffset)) // advance in block
			br.intermediate.Reset()
			for { // read sam records and append to answer until past the query region
				var curr Sam
				_, err = DecodeBam(br, &curr)
				if err == io.EOF {
					break
				}
				if curr.RName == chrom && curr.GetChromEnd() > int(start) && curr.GetChromStart() < int(end) {
					ans = append(ans, curr)
				}
				if (curr.RName == chrom && curr.GetChromStart() >= int(end)) || curr.RName != chrom {
					break
				}
			}
		}
	}
	return ans
}

// regionToBins returns all bins in which reads overlapping the input region may exist.
// Adapted from C code in SAM file specifications.
func regionToBins(beg, end int) []int {
	if beg == -1 && end == 0 {
		return []int{0, 0, 8, 72, 584, 4680}
	}
	var ans []int
	var k int
	end--
	for k = 1 + (beg >> 26); k <= 1+(end>>26); k++ {
		ans = append(ans, k)
	}
	for k = 9 + (beg >> 23); k <= 9+(end>>23); k++ {
		ans = append(ans, k)
	}
	for k = 73 + (beg >> 20); k <= 73+(end>>20); k++ {
		ans = append(ans, k)
	}
	for k = 585 + (beg >> 17); k <= 585+(end>>17); k++ {
		ans = append(ans, k)
	}
	for k = 4681 + (beg >> 14); k <= 4681+(end>>14); k++ {
		ans = append(ans, k)
	}
	return ans
}
