package interval

import (
	"github.com/vertgenlab/gonomics/axt"
	"github.com/vertgenlab/gonomics/bed"
	"github.com/vertgenlab/gonomics/chain"
	"github.com/vertgenlab/gonomics/sam"
	"github.com/vertgenlab/gonomics/vcf"
	"github.com/vertgenlab/gonomics/numbers"
	"path"
)

func LiftIntervalWithChain(c *chain.Chain, i Lift) (string, int, int) {
	return c.QName, chain.TPosToQPos(c, i.GetChromStart()), chain.TPosToQPos(c, i.GetChromEnd())
}

func GoReadToLiftChan(inputFile string) <-chan Lift {
	answer := make(chan Lift, 1000)
	go ReadToLiftChan(inputFile, answer)
	return answer
}

func ReadToLiftChan(inputFile string, send chan<- Lift) {
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
			send <- val
		}
	}
	/* These data types are theoretically compatible with Lift but not fully implemented.
	case ".axt":
		receive := axt.GoReadToChan(inputFile)
		for val := range receive {
			send <- val
		}

	case ".sam":
		receive, _ := sam.GoReadToChan(inputFile)
		for val := range receive {
			send <- val
		}
	case ".chain":
		receive, _ := chain.GoReadToChan(inputFile)
		for val := range receive {
			send <- val
		}
	*/
	close(send)
}

//matchOverlapLen returns the number of bases shared between two start and end points. Used in MatchProportion.
func matchOverlapLen(start1 int, end1 int, start2 int,  end2 int) int {
	return numbers.Max(0, numbers.Min(end1, end2)-numbers.Max(start1, start2))
}

//match proportion returns the proportion of bases in the target and query that can be lifted for a particular interval as a pair of floats (propT, propQ)
func MatchProportion(c *chain.Chain, i Interval) (float64, float64) {
	var match, dT, dQ int = 0, 0, 0
	var currPos int = c.TStart//starting with strand +/+ case for now.
	for j := 0; j < len(c.Alignment); j++ {
		match += matchOverlapLen(currPos, currPos + c.Alignment[j].Size, i.GetChromStart(), i.GetChromEnd())
		currPos += c.Alignment[j].Size
		dT += matchOverlapLen(currPos, currPos + c.Alignment[j].TBases, i.GetChromStart(), i.GetChromEnd())
		if matchOverlapLen(currPos, currPos + c.Alignment[j].TBases, i.GetChromStart(), i.GetChromEnd()) > 0 {
			dQ += c.Alignment[j].QBases//this handles a special case where both TBases and QBases are non-zero, as in competing non-aligned bases are present at the same location in target and query.
		}
		currPos += c.Alignment[j].TBases
	}
	return float64(match)/float64(match+dT), float64(match)/float64(match+dQ)
}

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
			send <- val
		}

	case ".sam":
		receive, _ := sam.GoReadToChan(inputFile)
		for val := range receive {
			send <- val
		}
	case ".chain":
		receive, _ := chain.GoReadToChan(inputFile)
		for val := range receive {
			send <- val
		}
	}
	close(send)
}
