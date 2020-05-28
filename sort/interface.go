package sort

import "github.com/vertgenlab/gonomics/fileio"

// Interface for generic merge sort of file types
type MergeSort interface {
	Len() int
	Swap(i, j int)
	Less(i, j int) bool
	Push(x interface{})
	Pop() interface{}
	Write(file string)
}

type MergeSortSingle interface {
	WriteToFileHandle(*fileio.EasyWriter)
	NextLine(*fileio.EasyReader) bool
}