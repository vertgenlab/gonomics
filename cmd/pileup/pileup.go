// Command Group: "SAM Tools"

package main

import (
	"flag"
	"fmt"
	"github.com/vertgenlab/gonomics/dna"
	"github.com/vertgenlab/gonomics/exception"
	"github.com/vertgenlab/gonomics/fileio"
	"github.com/vertgenlab/gonomics/sam"
	"log"
	"strings"
)

func usage() {
	fmt.Print(
		"pileup - Count bases from sequencing data\n\n" +
			"pileup prints to stdout by default, but can be redirected with the -o flag.\n" +
			"The output format is a tab-separated file with columns defined as:\n" +
			"Chr\tPos\t#A\t#C\t#G\t#T\t#N\t#DEL\tINS_array\n" +
			"INS_array is an array of observed insertions with format insSeq:readCount\n\n" +
			"Usage:\n" +
			"  pileup [options] in.bam\n\n" +
			"Options:\n")
	flag.PrintDefaults()
}

type Settings struct {
	minDp int
}

// pileFilters parses the Settings struct to make functions to filter the output piles.
func pileFilters(s Settings) []func(pile sam.Pile) bool {
	var filters []func(pile sam.Pile) bool

	if s.minDp > 0 {
		filters = append(filters, func(pile sam.Pile) bool {
			var count int
			for i := range pile.CountF {
				count += pile.CountF[i] + pile.CountR[i]
			}
			for _, i := range pile.InsCountF {
				count += i
			}
			for _, i := range pile.InsCountR {
				count += i
			}
			return count >= s.minDp
		})
	}

	return filters
}

func pileup(infile string, outfile string, settings Settings) {
	sBuilder := new(strings.Builder)
	reads, recycle, header := sam.GoReadToChanRecycle(infile, 1000)
	sendChan := make(chan sam.Sam) // no buffer to satisfy recycle contract

	// goroutine sends from the read channel to an unbuffered channel
	// then returns the previous struct to the recycle channel.
	go func(<-chan *sam.Sam, chan<- *sam.Sam, chan<- sam.Sam) {
		var prevRead *sam.Sam
		for read := range reads {
			sendChan <- *read
			if prevRead != nil {
				recycle <- prevRead
			}
			prevRead = read
		}
		close(sendChan)
	}(reads, recycle, sendChan)

	pileChan := sam.GoPileup(sendChan, header, false, nil, pileFilters(settings))

	var err error
	output := fileio.EasyCreate(outfile)
	_, err = fmt.Fprintln(output, "#Chr\tPos\tA\tC\tG\tT\tN\tDEL\tINS")
	exception.PanicOnErr(err)
	for pile := range pileChan {
		_, err = fmt.Fprintln(output, fmtOutput(pile, header, sBuilder))
		exception.PanicOnErr(err)
	}
	err = output.Close()
	exception.PanicOnErr(err)
}

// fmtOutput formats each bases pile for writing. Subject to change.
func fmtOutput(pile sam.Pile, header sam.Header, builder *strings.Builder) string {
	builder.Reset()
	builder.WriteString(fmt.Sprintf("%s\t%d\t%d\t%d\t%d\t%d\t%d\t%d", header.Chroms[pile.RefIdx].Name, pile.Pos,
		pile.CountF[dna.A]+pile.CountR[dna.A], pile.CountF[dna.C]+pile.CountR[dna.C], pile.CountF[dna.G]+pile.CountR[dna.G],
		pile.CountF[dna.T]+pile.CountR[dna.T], pile.CountF[dna.N]+pile.CountR[dna.N], pile.CountF[dna.Gap]+pile.CountR[dna.Gap]))
	for seq := range pile.InsCountF {
		builder.WriteString(fmt.Sprintf("\t%s:%d", seq, pile.InsCountF[seq]+pile.InsCountR[seq]))
	}
	var presentInForward bool
	for seq := range pile.InsCountR {
		_, presentInForward = pile.InsCountF[seq]
		if !presentInForward {
			builder.WriteString(fmt.Sprintf("\t%s:%d", seq, pile.InsCountR[seq]))
		}
	}
	return builder.String()
}

func main() {
	var output *string = flag.String("o", "stdout", "Output file")
	var minDp *int = flag.Int("minDP", 0, "Exclude positions with depth < minDP.")
	flag.Parse()
	flag.Usage = usage
	log.SetFlags(log.Ldate | log.Ltime)

	if len(flag.Args()) != 1 {
		usage()
		return
	}

	infile := flag.Arg(0)

	s := Settings{
		minDp: *minDp,
	}

	pileup(infile, *output, s)
}
