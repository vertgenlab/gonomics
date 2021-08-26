// Command Group: "FASTQ Tools"

package main

import (
	"flag"
	"fmt"
	"github.com/vertgenlab/gonomics/common"
	"github.com/vertgenlab/gonomics/exception"
	"github.com/vertgenlab/gonomics/fastq"
	"github.com/vertgenlab/gonomics/fileio"
	"github.com/vertgenlab/gonomics/numbers"
	"log"
	"math/rand"
	"sync"
)

type Settings struct {
	InFile        string
	OutFile       string
	R1InFile      string
	R2InFile      string
	R1OutFile     string
	R2OutFile     string
	PairedEnd     bool
	SubSet        float64
	RandSeed      bool
	SetSeed       int64
	MinSize       int
	MaxSize       int
	NamesList     string
	CollapseUmi   bool
	BarcodeLength int
	UmiLength     int
}

func fastqFilter(s Settings) {
	common.RngSeed(s.RandSeed, s.SetSeed)
	var r float64
	var err error
	var doneReading, FwdInMap, RevInMap, NameInMap bool
	var line string
	var currSc fastq.SingleCellPair
	namesMap := make(map[string]bool)

	if s.NamesList != "" {
		names := fileio.EasyOpen(s.NamesList)
		for line, doneReading = fileio.EasyNextRealLine(names); !doneReading; line, doneReading = fileio.EasyNextRealLine(names) {
			namesMap[line] = true
		}
		err = names.Close()
		exception.PanicOnErr(err)
	}

	if s.PairedEnd {
		umiMap := make(map[string]int)
		ReadCh := make(chan fastq.PairedEnd, 100000)
		WriteCh := make(chan fastq.PairedEnd, 100000)
		var wg sync.WaitGroup
		wg.Add(1)
		go fastq.PairedEndToChan(s.R1InFile, s.R2InFile, ReadCh)
		go fastq.WritingChan(s.R1OutFile, s.R2OutFile, WriteCh, &wg)
		for i := range ReadCh {
			if len(i.Fwd.Seq) < s.MinSize {
				continue
			}
			if len(i.Rev.Seq) < s.MinSize {
				continue
			}
			if len(i.Fwd.Seq) > s.MaxSize {
				continue
			}
			if len(i.Rev.Seq) > s.MaxSize {
				continue
			}
			if s.SubSet < 1 {
				r = rand.Float64()
				if r > s.SubSet {
					continue
				}
			}
			if s.NamesList != "" {
				if _, FwdInMap = namesMap[i.Fwd.Name]; !FwdInMap {
					if _, RevInMap = namesMap[i.Rev.Name]; !RevInMap {
						continue //continue if neither the forward or reverse read is in the map
					}
				}
			}
			if s.CollapseUmi {
				currSc = fastq.PairedEndToSingleCellPair(i, s.BarcodeLength, s.UmiLength)
				if fastq.UmiCollision(currSc, umiMap) {
					continue
				}
			}
			WriteCh <- i
		}
		close(WriteCh)
		wg.Wait()
	} else {
		f := fastq.GoReadToChan(s.InFile)
		out := fileio.EasyCreate(s.OutFile)

		for i := range f {
			r = rand.Float64()
			if r > s.SubSet {
				continue
			}
			if len(i.Seq) < s.MinSize {
				continue
			}
			if len(i.Seq) > s.MaxSize {
				continue
			}
			if s.NamesList != "" {
				if _, NameInMap = namesMap[i.Name]; !NameInMap {
					continue
				}
			}
			fastq.WriteToFileHandle(out, i)
		}
		err = out.Close()
		exception.PanicOnErr(err)
	}
}

func usage() {
	fmt.Print(
		"fastqFilter - Returns a filtered fastq based on option parameters.\n" +
			"Usage:\n" +
			"fastqFilter input.fastq output.fastq\n" +
			"OR\n" +
			"fastqFilter -pairedEnd R1.fastq R2.fastq out1.fastq out2.fastq\n" +
			"options:\n")
	flag.PrintDefaults()
}

func main() {
	var expectedNumArgs int = 2
	var pairedEnd *bool = flag.Bool("pairedEnd", false, "Paired end reads, use two input and output fastq files.")
	var subSet *float64 = flag.Float64("subSet", 1.0, "Proportion of reads to retain in output.")
	var randSeed *bool = flag.Bool("randSeed", false, "Uses a random seed for the RNG.")
	var setSeed *int64 = flag.Int64("setSeed", -1, "Use a specific seed for the RNG.")
	var minSize *int = flag.Int("minSize", 0, "Retain fastq reads above this size.")
	var maxSize *int = flag.Int("maxSize", numbers.MaxInt, "Retain fastq reads below this size.")
	var namesList *string = flag.String("namesList", "", "Specifies a file name for a line-delimited file of read names. Only reads with names matching the names in this file will be retained. For paired end reads, the read will be retained if either the forward or reverse read has a matching read name.")
	var collapseUmi *bool = flag.Bool("collapseUmi", false, "Removes UMI duplicates from single-cell format reads. R1: barcode and UMI. R2: mRNA sequence.")
	var barcodeLength *int = flag.Int("barcodeLength", 16, "Sets the length of the barcode for single-cell reads.")
	var umiLength *int = flag.Int("umiLength", 12, "Sets the UMI length for single-cell reads.")

	flag.Usage = usage
	log.SetFlags(log.Ldate | log.Ltime | log.Lshortfile)
	flag.Parse()

	if *subSet < 0 || *subSet > 1 {
		log.Fatalf(fmt.Sprintf("The subSet option must be between 0 and 1, received %v.", *subSet))
	}

	if *pairedEnd {
		expectedNumArgs = 4
	}

	if *collapseUmi && !*pairedEnd {
		log.Fatalf("To collapse UMIs from single-cell reads, select pairedEnd AND collapseUmi.")
	}

	if len(flag.Args()) != expectedNumArgs {
		flag.Usage()
		log.Fatalf("Error: expecting %d arguments, but got %d\n",
			expectedNumArgs, len(flag.Args()))
	}

	s := Settings{
		InFile:        "",
		OutFile:       "",
		R1InFile:      "",
		R2InFile:      "",
		R1OutFile:     "",
		R2OutFile:     "",
		PairedEnd:     *pairedEnd,
		SubSet:        *subSet,
		RandSeed:      *randSeed,
		SetSeed:       *setSeed,
		MinSize:       *minSize,
		MaxSize:       *maxSize,
		NamesList:     *namesList,
		CollapseUmi:   *collapseUmi,
		BarcodeLength: *barcodeLength,
		UmiLength:     *umiLength,
	}

	if *pairedEnd {
		s.R1InFile = flag.Arg(0)
		s.R2InFile = flag.Arg(1)
		s.R1OutFile = flag.Arg(2)
		s.R2OutFile = flag.Arg(3)
	} else {
		s.InFile = flag.Arg(0)
		s.OutFile = flag.Arg(1)
	}

	fastqFilter(s)
}
