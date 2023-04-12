// Command Group: "Sorting"

package main

import (
	"flag"
	"fmt"
	"log"
	"path"
	"strings"

	"github.com/vertgenlab/gonomics/axt"
	"github.com/vertgenlab/gonomics/bed"
	"github.com/vertgenlab/gonomics/dna"
	"github.com/vertgenlab/gonomics/exception"
	"github.com/vertgenlab/gonomics/fileio"
	"github.com/vertgenlab/gonomics/sam"
	"github.com/vertgenlab/gonomics/sort"
	"github.com/vertgenlab/gonomics/vcf"
)

func usage() {
	fmt.Print(
		"mergesort - Executes an external merge sort of the input file based on desired sort criteria. \n" +
			"\t The input file should have a proper file extension depending on the input file type.\n" +
			"\tDefault sort criteria is byGenomicCoordinates. Chromosome -> StartPos -> EndPos\n" +
			"Usage:\n" +
			" mergesort [options] input.filetype outputFile\n")
	flag.PrintDefaults()
}

func axtSort(infile, outfile string, numRecordsPerChunk int) {
	data, header := axt.GoReadToChan(infile)
	out := sort.GoExternalMergeSort(data, numRecordsPerChunk, func(a, b axt.Axt) bool {
		switch {
		case a.RName < b.RName:
			return true
		case a.RName > b.RName:
			return false
		case a.RStart < b.RStart:
			return true
		case a.RStart > b.RStart:
			return false
		default:
			return a.REnd < b.REnd
		}
	})

	o := fileio.EasyCreate(outfile)
	if len(header) != 0 {
		_, err := fmt.Fprintln(o, strings.Join(header, "\n"))
		exception.PanicOnErr(err)
	}
	var i int
	for r := range out {
		axt.WriteToFileHandle(o, r, i)
		i++
	}

	err := o.Close()
	exception.PanicOnErr(err)
}

func bedSort(infile, outfile string, numRecordsPerChunk int) {
	data := bed.GoReadToChan(infile)
	out := sort.GoExternalMergeSort(data, numRecordsPerChunk, func(a, b bed.Bed) bool {
		switch {
		case a.Chrom < b.Chrom:
			return true
		case a.Chrom > b.Chrom:
			return false
		case a.ChromStart < b.ChromStart:
			return true
		case a.ChromStart > b.ChromStart:
			return false
		default:
			return a.ChromEnd < b.ChromEnd
		}
	})

	o := fileio.EasyCreate(outfile)
	for r := range out {
		bed.WriteToFileHandle(o, r)
	}

	err := o.Close()
	exception.PanicOnErr(err)
}

func vcfSort(infile, outfile string, numRecordsPerChunk int) {
	data, header := vcf.GoReadToChan(infile)
	out := sort.GoExternalMergeSort(data, numRecordsPerChunk, func(a, b vcf.Vcf) bool {
		switch {
		case a.Chr < b.Chr:
			return true
		case a.Chr > b.Chr:
			return false
		default:
			return a.Pos < b.Pos
		}
	})

	o := fileio.EasyCreate(outfile)
	if len(header.Text) != 0 {
		_, err := fmt.Fprintln(o, strings.Join(header.Text, "\n"))
		exception.PanicOnErr(err)
	}
	for r := range out {
		vcf.WriteVcf(o, r)
	}

	err := o.Close()
	exception.PanicOnErr(err)
}

func samSort(infile, outfile string, numRecordsPerChunk int, sortCriteria string) {
	data, header := sam.GoReadToChan(infile)
	var out <-chan sam.Sam
	if sortCriteria == "singleCellBx" {
		out = sort.GoExternalMergeSort(data, numRecordsPerChunk, func(a, b sam.Sam) bool {
			iSingle := sam.ToSingleCellAlignment(a)
			jSingle := sam.ToSingleCellAlignment(b)
			if dna.BasesToString(iSingle.Bx) < dna.BasesToString(jSingle.Bx) {
				return true
			}
			return false
		})
	} else {
		out = sort.GoExternalMergeSort(data, numRecordsPerChunk, func(a, b sam.Sam) bool {
			switch {
			case a.RName < b.RName:
				return true
			case a.RName > b.RName:
				return false
			default:
				return a.Pos < b.Pos
			}
		})
	}

	o := fileio.EasyCreate(outfile)
	if len(header.Text) != 0 {
		_, err := fmt.Fprintln(o, strings.Join(header.Text, "\n"))
		exception.PanicOnErr(err)
	}
	for r := range out {
		sam.WriteToFileHandle(o, r)
	}

	err := o.Close()
	exception.PanicOnErr(err)
}

// TODO remove giraf pointers and uncomment
//func girafSort(infile, outfile string, numRecordsPerChunk int) {
//	data := giraf.GoReadToChan(infile)
//	out := sort.GoExternalMergeSort(data, numRecordsPerChunk, func(a, b *giraf.Giraf) bool {
//		// First sort criteria is node
//		if a.GetChrom() < b.GetChrom() {
//			return true
//		} else if a.GetChrom() == b.GetChrom() {
//			// If start nodes are equal then sort by start position
//			if a.GetChromStart() < b.GetChromStart() {
//				return true
//			} else if a.GetChromStart() == b.GetChromStart() {
//				// If start positions are equal then loop through nodes and see if one has priority
//				minPathLength := len(a.Path.Nodes)
//				if len(b.Path.Nodes) < minPathLength {
//					minPathLength = len(b.Path.Nodes)
//				}
//				for k := 0; k < minPathLength; k++ {
//					if a.Path.Nodes[k] < b.Path.Nodes[k] {
//						return true
//					}
//				}
//				// If all nodes match, sort based on longest path
//				if len(a.Path.Nodes) < len(b.Path.Nodes) {
//					return true
//				} else if len(a.Path.Nodes) == len(b.Path.Nodes) {
//					// If nodes are equal length, then sort based on the ending position
//					if a.GetChromEnd() < b.GetChromEnd() {
//						return true
//					}
//				}
//			}
//		}
//		return false
//	})
//
//	o := fileio.EasyCreate(outfile)
//	for r := range out {
//		giraf.WriteGirafToFileHandle(o, r)
//	}
//
//	err := o.Close()
//	exception.PanicOnErr(err)
//}

func mergeSort(filename string, outFile string, numRecordsPerChunk int, sortCriteria string) {
	// How the file is read is dependent on the file extension
	filetype := path.Ext(filename)

	if filetype == ".gz" {
		// If terminal extension is ".gz" then trim off the gz and get the next extension
		filetype = path.Ext(filename[0 : len(filename)-len(filetype)])
	}

	switch filetype {
	case ".axt":
		axtSort(filename, outFile, numRecordsPerChunk)
	case ".bed":
		bedSort(filename, outFile, numRecordsPerChunk)
	case ".vcf":
		vcfSort(filename, outFile, numRecordsPerChunk)
	case ".sam":
		samSort(filename, outFile, numRecordsPerChunk, sortCriteria)
	case ".giraf":
		// TODO enable after giraf pointers are removed
		log.Fatalln("ERROR: giraf sorting in currently disabled")
		//girafSort(filename, outFile, numRecordsPerChunk)
	default:
		if filetype == "" {
			log.Fatalln("ERROR: Input file must have proper file extension")
		} else {
			log.Fatalln("ERROR: Merge sort methods have not been implemented for file type:", filetype)
		}
	}
	return
}

func main() {
	expectedNumArgs := 2
	var numLinesPerChunk *int = flag.Int("tmpsize", 1000000, "The number of records to read into memory before writing to a tmp file.``")
	var singleCellBx *bool = flag.Bool("singleCellBx", false, "Sort single-cell sam records by barcode.")
	var sortCriteria string = "byGenomicCoordinates" //default the genomicCoordinates criteria.
	flag.Usage = usage
	log.SetFlags(log.Ldate | log.Ltime | log.Lshortfile)
	flag.Parse()

	if *singleCellBx {
		sortCriteria = "singleCellBx"
	}

	if len(flag.Args()) != expectedNumArgs {
		flag.Usage()
		log.Fatalf("ERROR: expecting %d arguments, but got %d\n", expectedNumArgs, len(flag.Args()))
	}

	inFile := flag.Arg(0)
	outFile := flag.Arg(1)

	mergeSort(inFile, outFile, *numLinesPerChunk, sortCriteria)
}
