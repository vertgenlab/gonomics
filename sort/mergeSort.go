package sort

import (
	"container/heap"
	"fmt"
	"github.com/vertgenlab/gonomics/axt"
	"github.com/vertgenlab/gonomics/bed"
	"github.com/vertgenlab/gonomics/common"
	"github.com/vertgenlab/gonomics/fileio"
	"github.com/vertgenlab/gonomics/giraf"
	"github.com/vertgenlab/gonomics/sam"
	"github.com/vertgenlab/gonomics/vcf"
	"log"
	"os"
	"path"
	"sort"
)

const (
	// Newly added filetypes should be added to the following valid types list
	validFileTypes     = "axt, bed, vcf, sam, giraf"
	maxTmpFilesAllowed = 1000
)

// New case must be added to readChunk for each new filetype
// This function declare a MergeSort and MergeSortSingle with the
// underlying type corresponding to the filetype input
func chooseDataType(filetype string, sortType string) (MergeSort, MergeSortSingle) {
	var answer MergeSort
	var single MergeSortSingle
	switch filetype {
	case ".axt":
		switch sortType {
		case "byGenomicCoordinates":
			answer = new(axt.ByGenomicCoordinates)
		default:
			log.Fatalf("ERROR: Sorting %s has not been implemented for %s", sortType, filetype)
		}
		single = new(axt.Axt)

	case ".bed":
		switch sortType {
		case "byGenomicCoordinates":
			answer = new(bed.ByGenomicCoordinates)
		default:
			log.Fatalf("ERROR: Sorting %s has not been implemented for %s", sortType, filetype)
		}
		single = new(bed.Bed)

	case ".vcf":
		switch sortType {
		case "byGenomicCoordinates":
			answer = new(vcf.ByGenomicCoordinates)
		default:
			log.Fatalf("ERROR: Sorting %s has not been implemented for %s", sortType, filetype)
		}
		single = new(vcf.Vcf)

	case ".sam":
		switch sortType {
		case "byGenomicCoordinates":
			answer = new(sam.ByGenomicCoordinates)
		default:
			log.Fatalf("ERROR: Sorting %s has not been implemented for %s", sortType, filetype)
		}
		single = new(sam.Sam)

	case ".giraf":
		switch sortType {
		case "byGenomicCoordinates":
			answer = new(giraf.ByTopologicalNodeOrder)
		default:
			log.Fatalf("ERROR: Sorting %s has not been implemented for %s", sortType, filetype)
		}
		single = new(giraf.Giraf)

	default:
		log.Println("Valid file types include:", validFileTypes)
		if filetype == "" {
			log.Fatalln("ERROR: Input file must have proper file extension")
		} else {
			log.Fatalln("ERROR: Merge sort methods have not been implemented for file type:", filetype)
		}
	}
	return answer, single
}

func ExternalMergeSort(filename string, linesPerChunk int, tmpFilePrefix string, outFilename string, sortType string) {
	// Default sort order is by genomic coordinates
	if sortType == "" {
		sortType = "byGenomicCoordinates"
	}

	// How the file is read is dependent on the file extension
	filetype := path.Ext(filename)

	if filetype == ".gz" {
		// If terminal extension is ".gz" then trim off the gz and get the next extension
		filetype = path.Ext(filename[0 : len(filename)-len(filetype)])
	}

	file := fileio.EasyOpen(filename)
	var done bool
	var tmpFiles []string
	var tmpFilename string
	var tmpFileId int = 0
	var currChunk MergeSort

	// Read files into sorted chunk until reached EOF
	for currChunk, done = readChunk(file, linesPerChunk, filetype, sortType); currChunk.Len() != 0; currChunk, done = readChunk(file, linesPerChunk, filetype, sortType) {
		if tmpFileId == maxTmpFilesAllowed {
			log.Fatalln("ERROR: Exceeded maximum number of tmp files, increase -chunkSize")
		}

		sort.Sort(currChunk) // Sort the incoming chunk to be written

		// Handle naming of tmp file and keep a record of it
		tmpFilename = fmt.Sprintf("%s_%d", tmpFilePrefix, tmpFileId)
		tmpFiles = append(tmpFiles, tmpFilename)
		tmpFileId++

		currChunk.Write(tmpFilename) // Write the sorted chunk to a tmp file

		if done {
			break
		}
	}

	mergeChunks(tmpFiles, outFilename, filetype, sortType) // Begin merge of tmp files
}

func readChunk(file *fileio.EasyReader, lines int, filetype string, sortType string) (MergeSort, bool) {
	var done bool = false
	// Initialize an empty chunk and an empty single value with the proper underlying type for later functions
	chunk, curr := chooseDataType(filetype, sortType)

	// In the following function the curr variable acts as a static pointer to a changing value that is
	// updated with each call of curr.NextRealLine. To actually get each value in a slice, we need
	// to copy the underlying value in the pointer to an interface (valueToAdd) then assert the type of
	// that copied value and push it onto the heap

	for i := 0; i < lines; i++ {
		done = curr.NextRealRecord(file) // Get the next value and store it in the curr receiver
		if done {
			done = true
			break
		} else if curr != nil {
			valueToAdd := curr.Copy()                  // Copy the curr value to a new memory address
			chunk.Push((valueToAdd).(MergeSortSingle)) // Push the copied pointer to the slice
		}
	}
	return chunk, done
}

func mergeChunks(tmpFiles []string, outFilename string, filetype string, sortType string) {
	fileReaders := make([]*fileio.EasyReader, len(tmpFiles))
	var done bool

	// This memory address map will store the memory address of each of the input variables
	// and key them to the origin file reader. The memory addresses for each file is kept static
	// throughout the heap loop below and therefore can be used as a stable key to determine
	// the file of origin
	memoryAddressMap := make(map[interface{}]*fileio.EasyReader)

	priorityQueue, curr := chooseDataType(filetype, sortType)

	for i := 0; i < len(tmpFiles); i++ {
		fileReaders[i] = fileio.EasyOpen(tmpFiles[i])
		done = curr.NextRealRecord(fileReaders[i])

		valueToAdd := curr.Copy()                   // Copy the curr value to a new memory address
		writeVal := (valueToAdd).(MergeSortSingle)  // Assert the type of the empty interface valueToAdd
		priorityQueue.Push(writeVal)                // Push the copied pointer to the heap
		memoryAddressMap[writeVal] = fileReaders[i] // Key the memory address to the origin file reader
	}

	// Initialize heap
	heap.Init(priorityQueue)

	var currVal MergeSortSingle
	var currFile *fileio.EasyReader
	outFile := fileio.EasyCreate(outFilename)
	defer outFile.Close()
	// Pop heap until all tmp files are exhausted
	for priorityQueue.Len() > 0 {
		done = false
		// Get minimum value from the heap, recall that the memory address of this value
		// has been keyed in the memoryAddressMap to the origin file reader
		currVal = heap.Pop(priorityQueue).(MergeSortSingle)
		currFile = memoryAddressMap[currVal]    // Get the origin file reader
		currVal.WriteToFileHandle(outFile)      // Write the lowest value to the outfile
		done = currVal.NextRealRecord(currFile) // Get the next value from the origin file

		if done { // if done then close and remove the tmp file
			currFile.Close()
			err := os.Remove(currFile.File.Name())
			common.ExitIfError(err)
		} else {
			heap.Push(priorityQueue, currVal) // push the new value onto the heap
		}
	}
}
