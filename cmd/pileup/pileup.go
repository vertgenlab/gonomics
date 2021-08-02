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
			"Usage:\n" +
			"  pileup [options] in.bam\n\n" +
			"Options:\n")
	flag.PrintDefaults()
}

type Settings struct {
	minDp int
}

// pileFilters parses the Settings struct to make functions to filter the output piles.
func pileFilters(s Settings) []func(p sam.Pile) bool {
	var filters []func(p sam.Pile) bool

	if s.minDp > 0 {
		filters = append(filters, func(p sam.Pile) bool {
			var count int
			for i := range p.Count {
				count += p.Count[i]
			}
			for _, i := range p.InsCount {
				count += i
			}
			return count >= s.minDp
		})
	}

	return filters
}

func pileup(infile string, outfile string, s Settings) {
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

	pileChan := sam.GoPileup(sendChan, header, false, nil, pileFilters(s))

	var err error
	output := fileio.EasyCreate(outfile)
	_, err = fmt.Fprintln(output, "Chr\tPos\tA\tC\tG\tT\tN\tDEL\tINS")
	exception.PanicOnErr(err)
	for pile := range pileChan {
		_, err = fmt.Fprintln(output, fmtOutput(pile, header))
		exception.PanicOnErr(err)
	}
	err = output.Close()
	exception.PanicOnErr(err)
}

// fmtOutput formats each bases pile for writing. Subject to change.
func fmtOutput(pile sam.Pile, h sam.Header) string {
	s.Reset()
	s.WriteString(fmt.Sprintf("%s\t%d\t%d\t%d\t%d\t%d\t%d\t%d", h.Chroms[pile.RefIdx].Name, pile.Pos,
		pile.Count[dna.A], pile.Count[dna.C], pile.Count[dna.G], pile.Count[dna.T], pile.Count[dna.N], pile.Count[dna.Gap]))
	for seq, count := range pile.InsCount {
		s.WriteString(fmt.Sprintf("\t%s:%d", seq, count))
	}
	return s.String()
}

var s *strings.Builder

func init() {
	s = new(strings.Builder)
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
