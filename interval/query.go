package interval

import (
	"github.com/vertgenlab/gonomics/axt"
	"github.com/vertgenlab/gonomics/bed"
	"github.com/vertgenlab/gonomics/chain"
	"github.com/vertgenlab/gonomics/sam"
	"github.com/vertgenlab/gonomics/vcf"
	"path"
)

func GoReadToChan(inputFile string) <-chan Interval {
	answer := make(chan Interval, 1000)
	go ReadToChan(inputFile, answer)
	return answer
}

func ReadToChan(inputFile string, send chan<- Interval) {
	// How the file is read is dependent on the file extension
	filetype := path.Ext(inputFile)

	if filetype == ".gz" {
		// If terminal extension is ".gz" then trim off the gz and get the next extension
		filetype = path.Ext(inputFile[0 : len(inputFile)-len(filetype)])
	}

	switch filetype {
	case ".bed":
		receive := bed.GoReadToChan(inputFile)
		for val := range receive {
			send <- val
		}

	case ".axt":
		receive := axt.GoReadToChan(inputFile)
		for val := range receive {
			send <- val
		}

	case ".vcf":
		receive, _ := vcf.GoReadToChan(inputFile)
		for val := range receive {
			send <- &val
		}

	case ".sam":
		receive, _ := sam.GoReadToChan(inputFile)
		for val := range receive {
			send <- &val
		}
	case ".chain":
		receive, _ := chain.GoReadToChan(inputFile)
		for val := range receive {
			send <- val
		}
	}
	close(send)
}
