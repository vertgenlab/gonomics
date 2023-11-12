// Command Group: "Sorting"

// Executes an external merge sort of the input file based on desired sort criteria
package main

import (
	"flag"
	"fmt"
	"github.com/vertgenlab/gonomics/fastq"
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
			"\tThe input file should have a proper file extension depending on the input file type.\n" +
			"\tDefault sort criteria is byGenomicCoordinates. Chromosome -> StartPos -> EndPos\n" +
			"\tExcept in the case of fastq files, where the file is sorted by read name.\n" +
			"\tCurrent accepted file types: BED, VCF, SAM, AXT, FASTQ\n" +
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

func samSortSpecial(inFile, outFile string, numRecordsPerChunk int) {
	var o <-chan sam.Sam
	//var fieldsA, fieldsB []string
	reads, head := sam.GoReadToChan(inFile)
	out := fileio.EasyCreate(outFile)
	sam.WriteHeaderToFileHandle(out, head)
	o = sort.GoExternalMergeSort(reads, numRecordsPerChunk, func(a, b sam.Sam) bool {
		return a.QName < b.QName
	})
	for i := range o {
		sam.WriteToFileHandle(out, i)
	}
	err := out.Close()
	exception.PanicOnErr(err)
}

func samSort(infile, outfile string, numRecordsPerChunk int, sortCriteria string) {
	data, header := sam.GoReadToChan(infile)
	var out <-chan sam.Sam
	if sortCriteria == "singleCellBx" {
		out = sort.GoExternalMergeSort(data, numRecordsPerChunk, func(a, b sam.Sam) bool {
			iSingle := sam.ToSingleCellAlignment(a)
			jSingle := sam.ToSingleCellAlignment(b)
			return dna.BasesToString(iSingle.Bx) < dna.BasesToString(jSingle.Bx)
		})
	} else if sortCriteria == "readName" {
		out = sort.GoExternalMergeSort(data, numRecordsPerChunk, func(a, b sam.Sam) bool {
			return a.QName < b.QName
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

func fastqSort(inFile string, outFile string, numRecordsPerChunk int) {
	//detect if PE mode
	peIn := strings.Split(inFile, ",")
	peOut := strings.Split(outFile, ",")

	for file := range peIn {
		var outChan <-chan fastq.Fastq

		out := fileio.EasyCreate(peOut[file])
		data := fastq.GoReadToChan(peIn[file])

		outChan = sort.GoExternalMergeSort(data, numRecordsPerChunk, func(a, b fastq.Fastq) bool {
			return a.Name < b.Name
		})

		for fq := range outChan {
			fastq.WriteToFileHandle(out, fq)
		}

		err := out.Close()
		exception.PanicOnErr(err)
	}
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
	case ".fastq":
		fastqSort(filename, outFile, numRecordsPerChunk)
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
}

func main() {
	var numLinesPerChunk *int = flag.Int("tmpsize", 1000000, "The number of records to read into memory before writing to a tmp file.``")
	var singleCellBx *bool = flag.Bool("singleCellBx", false, "Sort single-cell sam records by barcode.")
	var sortCriteria string = "byGenomicCoordinates" // default the genomicCoordinates criteria.
	var fastqPE *bool = flag.Bool("fastqPE", false, "Sort paired end fastq files. Each file will be sorted separately")
	var readName *bool = flag.Bool("readName", false, "Sort sam records by read name.")
	flag.Usage = usage
	log.SetFlags(log.Ldate | log.Ltime | log.Lshortfile)
	flag.Parse()

	var expectedNumArgs int = 2
	if *fastqPE {
		expectedNumArgs = 4
	}

	if *singleCellBx {
		sortCriteria = "singleCellBx"
	}
	if *readName {
		sortCriteria = "readName"
	}

	if len(flag.Args()) != expectedNumArgs {
		flag.Usage()
		log.Fatalf("ERROR: expecting %d arguments, but got %d\n", expectedNumArgs, len(flag.Args()))
	}

	var inFile, outFile string
	if *fastqPE {
		inFile = fmt.Sprintf("%s,%s", flag.Arg(0), flag.Arg(1))
		outFile = fmt.Sprintf("%s,%s", flag.Arg(2), flag.Arg(3))
	} else {
		inFile = flag.Arg(0)
		outFile = flag.Arg(1)
	}

	mergeSort(inFile, outFile, *numLinesPerChunk, sortCriteria)
}
