package sort

import (
	"container/heap"
	"fmt"
	"github.com/vertgenlab/gonomics/axt"
	"github.com/vertgenlab/gonomics/bed"
	"github.com/vertgenlab/gonomics/common"
	"github.com/vertgenlab/gonomics/fileio"
	"github.com/vertgenlab/gonomics/vcf"
	"log"
	"os"
	"path"
	"sort"
)

const (
	// Newly added filetypes should be added to the following valid types list
	validFileTypes = "axt, bed, vcf"
	maxTmpFilesAllowed = 1000
)

func ExternalMergeSort(filename string, linesPerChunk int, tmpFilePrefix string, outFilename string) {

	// How the file is read is dependent on the file extension
	filetype := path.Ext(filename)

	file := fileio.EasyOpen(filename)
	var done bool
	var tmpFiles []string
	var tmpFilename string
	var tmpFileId int = 0
	var currChunk MergeSort

	// Read files into sorted chunk until reached EOF
	for currChunk, done = readChunk(file, linesPerChunk, filetype); !done; currChunk, done = readChunk(file, linesPerChunk, filetype) {
		if currChunk.Len() == 0 {
			continue
		}

		// Sort the incoming chunk to be written
		sort.Sort(currChunk)

		// Handle naming of tmp file and keep a record of it
		tmpFilename = fmt.Sprintf("%s_%d", tmpFilePrefix, tmpFileId)
		tmpFiles = append(tmpFiles, tmpFilename)
		tmpFileId++

		// Write the sorted chunk to a tmp file
		currChunk.Write(tmpFilename)
	}

	// Begin merge of tmp files
	mergeChunks(tmpFiles, outFilename, filetype)
}

// New case must be added to readChunk for each new filetype
func readChunk(file *fileio.EasyReader, lines int, filetype string) (MergeSort, bool) {
	var done bool = false
	switch filetype {

	// Functions to read in each filetype line by line
	// New case must be added for each new file type to be sorted

	// For the most part these can be copy-paste so long as the
	// NextLine function has been implemented for the desired data type

	case "axt":
		var curr *axt.Axt
		chunk := make([]*axt.Axt, 0, lines)
		for i := 0; i < lines; i++ {
			curr, done = axt.NextAxt(file)
			if done {
				done = true
				break
			}
			chunk = append(chunk, curr)
		}
		answer := axt.ByGenomicCoordinates(chunk)
		return &answer, done

	case "bed":
		var curr *bed.Bed
		chunk := make([]*bed.Bed, 0, lines)
		for i := 0; i < lines; i++ {
			curr, done = bed.NextBed(file)
			if done {
				done = true
				break
			}
			chunk = append(chunk, curr)
		}
		answer := bed.ByGenomicCoordinates(chunk)
		return &answer, done

	case "vcf":
		var curr *vcf.Vcf
		chunk := make([]*vcf.Vcf, 0, lines)
		for i := 0; i < lines; i++ {
			curr, done = vcf.NextVcf(file)
			if done {
				done = true
				break
			}
			chunk = append(chunk, curr)
		}
		answer := vcf.ByGenomicCoordinates(chunk)
		return &answer, done

	default:
		log.Println("Valid file types include:", validFileTypes)
		if filetype == "" {
			log.Fatalln("ERROR: Input file must have proper file extension")
		} else {
			log.Fatalln("ERROR: Merge sort methods have not been implemented for file type:", filetype)
		}
	}
	// Should never get to this point
	return nil, true
}

// New case must be added to mergeChunk for each new filetype
func mergeChunks(tmpFiles []string, outFilename string, filetype string) {
	fileReaders := make([]*fileio.EasyReader, len(tmpFiles))
	var priorityQueue MergeSort
	var done bool

	// Stores the memory address of variable from each file and retrieves the file it came from
	memoryAddressMap := make(map[interface{}]*fileio.EasyReader)

	// Must add new filetype to list AND a new case
	// TODO: better way to handle this???
	var Axt axt.ByGenomicCoordinates
	var Bed bed.ByGenomicCoordinates
	var Vcf vcf.ByGenomicCoordinates

	for i := 0; i < len(tmpFiles); i++ {
		fileReaders[i] = fileio.EasyOpen(tmpFiles[i])

		switch filetype {

		case "axt":
			var curr *axt.Axt
			curr, done = axt.NextAxt(fileReaders[i])
			memoryAddressMap[curr] = fileReaders[i]
			if !done{Axt = append(Axt, curr)}

		case "bed":
			var curr *bed.Bed
			curr, done = bed.NextBed(fileReaders[i])
			memoryAddressMap[curr] = fileReaders[i]
			if !done{Bed = append(Bed, curr)}

		case "vcf":
			var curr *vcf.Vcf
			curr, done = vcf.NextVcf(fileReaders[i])
			memoryAddressMap[curr] = fileReaders[i]
			if !done{Vcf = append(Vcf, curr)}
		}
		if done {
			fileReaders[i].Close()
			err := os.Remove(fileReaders[i].File.Name())
			common.ExitIfError(err)
			log.Fatalln("ERROR: Could not find file, or file empty", tmpFiles[i])
		}
	}


	switch filetype {
	case "axt":
		priorityQueue = &Axt

	case "bed":
		priorityQueue = &Bed

	case "vcf":
		priorityQueue = &Vcf
	}

	// Initialize heap
	heap.Init(priorityQueue)

	// Pop heap until all tmp files are exhausted
	var curr MergeSortSingle
	var currFile *fileio.EasyReader
	outFile := fileio.EasyCreate(outFilename)
	defer outFile.Close()
	for priorityQueue.Len() < 0 {
		done = false
		curr = heap.Pop(priorityQueue).(MergeSortSingle)
		currFile = memoryAddressMap[curr]
		curr.WriteToFileHandle(outFile)
		done = curr.NextLine(currFile)

		if done {
			currFile.Close()
			err := os.Remove(currFile.File.Name())
			common.ExitIfError(err)
		} else {
			heap.Push(priorityQueue, curr)
		}
	}
}