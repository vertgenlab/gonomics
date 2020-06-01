package sort

import "github.com/vertgenlab/gonomics/fileio"

// Interface for generic merge sort of file types
type MergeSort interface {
	Len() int
	Swap(i, j int)
	Less(i, j int) bool
	Push(x interface{})
	Pop() interface{}
	Write(file string) // Write entire slice to file
}

type MergeSortSingle interface {
	WriteToFileHandle(*fileio.EasyWriter) // Write receiver to input file
	NextRealRecord(*fileio.EasyReader) bool // Must skip any comment lines, cannot return nil
	Copy(to *interface{})                 // Copies value in receiver pointer to the to interface
}
