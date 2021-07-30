// Command Group: "SAM Tools"

package main

import (
	"flag"
	"fmt"
	"github.com/vertgenlab/gonomics/sam"
	"log"
)

func usage() {
	fmt.Print(
		"pileup - Count bases from sequencing data\n\n" +
			"Usage:\n" +
			"  pileup [options] in.bam\n\n" +
			"Options:\n")
	flag.PrintDefaults()
}

func pileup(infile string) {
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

	pileChan := sam.GoPileup(sendChan, header, nil)

	for pile := range pileChan {
		fmt.Println(pile)
	}
}

func main() {
	var output *string = flag.String("o", "stdout", "output file")
	flag.Parse()
	flag.Usage = usage
	log.SetFlags(log.Ldate | log.Ltime)

	if len(flag.Args()) != 1 {
		usage()
		return
	}

	pileup(*output)
}
