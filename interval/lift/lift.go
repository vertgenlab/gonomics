package lift

import (
	"github.com/vertgenlab/gonomics/bed"
	"github.com/vertgenlab/gonomics/chain"
	"github.com/vertgenlab/gonomics/interval"
	"github.com/vertgenlab/gonomics/vcf"
	"io"
	"log"
	"path"
)

//Lift is an interface for genomic regions. Unlike interval, Lifts can be edited in place.
type Lift interface {
	GetChrom() string
	GetChromStart() int
	GetChromEnd() int
	WriteToFileHandle(io.Writer)
	UpdateCoord(string, int, int) interface{}
}

//GoRead reads Lift interfaces from an input file to a slice. Uses Chans under the hood.
func GoRead(inputFile string) []Lift {
	var answer []Lift
	ch := GoReadToChan(inputFile)
	for i := range ch {
		answer = append(answer, i)
	}
	return answer
}

//GoReadToChan reads Lift interfaces to a channel from an input file.
func GoReadToChan(inputFile string) <-chan Lift {
	answer := make(chan Lift, 1000)
	go ReadToChan(inputFile, answer)
	return answer
}

//ReadToChan reads from a file to send Lift interfaces to a chan<- Lift.
func ReadToChan(inputFile string, send chan<- Lift) {
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

	case ".vcf":
		receive, _ := vcf.GoReadToChan(inputFile)
		for val := range receive {
			curr := val
			send <- curr
		}
	default:
		log.Fatalf("Filetype does not satisfy Lift interface: %v.", inputFile)
	}
	close(send)
}

//LiftCoordinatesWithChain returns the update coordinates for an input Lift interface based on an input chain.
func LiftCoordinatesWithChain(c *chain.Chain, i Lift) (string, int, int) {
	var newStart, newEnd int
	newStart = chain.TPosToQPos(c, i.GetChromStart())
	/* The minus one/plus one handles when a region ends
	with a structural variant and ensures correct placement of the end in the new assembly. */
	newEnd = chain.TPosToQPos(c, i.GetChromEnd()-1) + 1
	if !c.QStrand { //valid bed formats must have start < end. So this correction is made for intervals lifted to the negative strand.
		newStart, newEnd = newEnd, newStart
		newStart += 1 //these lines are corrections for the open/closed interval start and end.
		newEnd += 1
	}
	return c.QName, newStart, newEnd
}

//IntervalSliceToLift casts a slice of intervals to a slice of Lift interfaces.
func IntervalSliceToLift(b []interval.Interval) []Lift {
	var answer []Lift = make([]Lift, len(b))
	for i := range b {
		answer[i] = b[i].(Lift)
	}
	return answer
}
