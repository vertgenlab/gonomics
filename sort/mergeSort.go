package sort

import (
	"container/heap"
	"encoding/gob"
	"github.com/vertgenlab/gonomics/exception"
	"io"
	"os"
	"sort"
)

type heapBuf []any

// lessFunc is used to define the Less method for heapBuf.
// lessFunc is set by the mergeSort function which inputs
// a less function for the input generic data type.
// It is not ideal to use a mutable global variable here,
// but I was not able to find a cleaner solution.
var lessFunc func(a, b any) bool

func (hb heapBuf) Len() int {
	return len(hb)
}

func (hb heapBuf) Swap(i, j int) {
	hb[i], hb[j] = hb[j], hb[i]
}

func (hb heapBuf) Less(i, j int) bool {
	return lessFunc(hb[i], hb[j])
}

func (hb *heapBuf) Push(x any) {
	*hb = append(*hb, x)
}

func (hb *heapBuf) Pop() any {
	val := (*hb)[len(*hb)-1]
	*hb = (*hb)[:len(*hb)-1]
	return val
}

// GoExternalMergeSort inputs generic data records from a channel and a function to compare multiple data records
// and outputs a channel of sorted data records. This function sorts externally via merge sort and therefore is
// appropriate for very large amounts of data.
func GoExternalMergeSort[E any](data <-chan E, recordsPerTmpFile int, less func(a, b E) bool) <-chan E {
	out := make(chan E, 1000)
	go mergeSort(data, less, out, recordsPerTmpFile)
	return out
}

func mergeSort[E any](data <-chan E, less func(a, b E) bool, out chan<- E, recordsPerTmpFile int) {
	var err error

	// write records from data channel into multiple temporary files (chunks)
	// that are each independently sorted using the input less function
	var tmpFiles []*os.File
	tmpFiles = chunkData(data, recordsPerTmpFile, less)

	// a bit of hacking to make types match so that we can use a generic method
	lessFunc = func(a, b any) bool {
		return less(*a.(*E), *b.(*E))
	}

	// This map will store the memory address of each index in a buffer slice (one per file)
	// and key them to the decoder for the origin file. The memory addresses for the buffer
	// is kept static throughout the heap loop below and therefore can be used as a stable key
	// to determine the file of origin
	memoryAddressMap := make(map[*E]*gob.Decoder)

	curr := make([]*E, len(tmpFiles)) // this functions as a buffer for hb with static memory addresses
	for i := range curr {
		curr[i] = new(E) // allocate mem for curr pointers
	}
	hb := make(heapBuf, len(tmpFiles))
	decoders := make([]*gob.Decoder, len(tmpFiles))

	for i := range hb {
		tmpFiles[i], err = os.Open(tmpFiles[i].Name())
		exception.PanicOnErr(err)
		decoders[i] = gob.NewDecoder(tmpFiles[i])
		err = decoders[i].Decode(curr[i])
		exception.PanicOnErr(err)
		hb[i] = curr[i]
		memoryAddressMap[curr[i]] = decoders[i]
	}

	// Initialize heap
	pq := &hb
	heap.Init(pq)

	// Pop heap until all tmp files are exhausted
	var currVal *E
	var empty E
	for hb.Len() > 0 {
		// Get minimum value from the heap, recall that the memory address of this value
		// has been keyed in the memoryAddressMap to the origin file reader
		currVal = heap.Pop(pq).(*E)
		out <- *currVal // send popped value
		// reset E to dealloc pointer refs in E so gob reallocs instead of overwrites.
		// If this is not done then Decode will overwrite references (e.g. slices) that
		// may be present in currVal leading to pointer errors.
		*currVal = empty
		err = memoryAddressMap[currVal].Decode(currVal) // read in a new value from the origin decoder
		if err == io.EOF {
			continue
		}
		exception.PanicOnErr(err)
		heap.Push(pq, currVal) // push the new value onto the heap
	}

	close(out)
}

// chunkData writes records from data to temporary files each containing recordsPerChunk number
// of records with each file independently sorted. The return is a slice of closed files that can
// used to retrieve the name for reopening.
func chunkData[E any](data <-chan E, recordsPerChunk int, less func(a, b E) bool) []*os.File {
	currChunk := make([]E, 0, recordsPerChunk)
	var tmpFiles []*os.File

	for v := range data {
		currChunk = append(currChunk, v)
		if len(currChunk) == recordsPerChunk {
			sort.Slice(currChunk, func(i, j int) bool {
				return less(currChunk[i], currChunk[j])
			})
			tmpFiles = append(tmpFiles, writeTmpFile(currChunk))
			currChunk = currChunk[:0]
		}
	}

	if len(currChunk) != 0 {
		sort.Slice(currChunk, func(i, j int) bool {
			return less(currChunk[i], currChunk[j])
		})
		tmpFiles = append(tmpFiles, writeTmpFile(currChunk))
		currChunk = currChunk[:0]
	}

	return tmpFiles
}

// writeTmpFile writes arbitrary type E to a file using the gob package for
// generating self-described binary that can be parsed into a struct by the
// gob decoder.
func writeTmpFile[E any](chunk []E) *os.File {
	file, err := os.CreateTemp("", "sort_chunk_")
	defer file.Close()
	exception.PanicOnErr(err)

	encoder := gob.NewEncoder(file)
	for i := range chunk {
		err = encoder.Encode(chunk[i]) // encode each data record
		exception.PanicOnErr(err)
	}
	return file
}
