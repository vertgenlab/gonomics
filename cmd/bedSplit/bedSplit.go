package main

import (
	"flag"
	"fmt"
	"github.com/vertgenlab/gonomics/bed"
	"github.com/vertgenlab/gonomics/exception"
	"github.com/vertgenlab/gonomics/fileio"
	"log"
	"os"
)

type Settings struct {
	Mode    string
	InFile  string
	OutDir  string
	GzipOut bool
}

func bedSplit(s Settings) {
	var err error
	var foundInMap bool
	if _, err = os.Stat(s.OutDir); os.IsNotExist(err) {
		err = os.Mkdir(s.OutDir, 0755)
		exception.PanicOnErr(err)
	}

	records := bed.GoReadToChan(s.InFile)
	switch s.Mode {
	case "byName":
		var seenNames = make(map[string]*fileio.EasyWriter, 0)
		for v := range records {
			if _, foundInMap = seenNames[v.Name]; !foundInMap {
				if s.GzipOut {
					seenNames[v.Name] = fileio.EasyCreate(fmt.Sprintf("%s/%s.bed.gz", s.OutDir, v.Name))
				} else {
					seenNames[v.Name] = fileio.EasyCreate(fmt.Sprintf("%s/%s.bed", s.OutDir, v.Name))
				}

			}
			bed.WriteBed(seenNames[v.Name], v)
		}
		for currFile := range seenNames {
			err = seenNames[currFile].Close()
			exception.PanicOnErr(err)
		}
	case "byChrom":
		var seenChroms = make(map[string]*fileio.EasyWriter, 0)
		for v := range records {
			if _, foundInMap = seenChroms[v.Chrom]; !foundInMap {
				if s.GzipOut {
					seenChroms[v.Chrom] = fileio.EasyCreate(fmt.Sprintf("%s/%s.bed.gz", s.OutDir, v.Chrom))
				} else {
					seenChroms[v.Chrom] = fileio.EasyCreate(fmt.Sprintf("%s/%s.bed", s.OutDir, v.Chrom))
				}
			}
			bed.WriteBed(seenChroms[v.Chrom], v)
		}
		for currFile := range seenChroms {
			err = seenChroms[currFile].Close()
			exception.PanicOnErr(err)
		}
	default:
		log.Fatalf("Error: unrecognized mode. Mode may be either 'byChrom' or 'byName'.\n")
	}
}

func usage() {
	fmt.Print(
		"bedSplit - partition an input bed file.\n" +
			"Usage:\n" +
			"bedSplit mode input.bed outDir\n" +
			"Mode may be one of the following: byName, byChrom\n" +
			"options:\n")
}

func main() {
	var expectedNumArgs int = 3
	var gzipOut *bool = flag.Bool("gzipOut", false, "Compress the output bed files in gzip format.")
	flag.Usage = usage
	log.SetFlags(log.Ldate | log.Ltime | log.Lshortfile)
	flag.Parse()

	if len(flag.Args()) != expectedNumArgs {
		flag.Usage()
		log.Fatalf("Error: expecting %d arguments, but got %d\n", expectedNumArgs, len(flag.Args()))
	}

	mode := flag.Arg(0)
	inFile := flag.Arg(1)
	outDir := flag.Arg(2)

	s := Settings{
		Mode:    mode,
		InFile:  inFile,
		OutDir:  outDir,
		GzipOut: *gzipOut,
	}

	bedSplit(s)
}
