// Package lift provides an interface and utilities for working with genomic regions and offer flexibility to modify their coordinates.
package lift

import (
	"github.com/vertgenlab/gonomics/axt"
	"github.com/vertgenlab/gonomics/dna"
	"io"
	"log"
	"path"

	"github.com/vertgenlab/gonomics/bed"
	"github.com/vertgenlab/gonomics/chain"
	"github.com/vertgenlab/gonomics/interval"
	"github.com/vertgenlab/gonomics/numbers"
	"github.com/vertgenlab/gonomics/vcf"
)

// Lift is an interface for genomic regions. Unlike interval, Lifts can be edited in place.
type Lift interface {
	GetChrom() string
	GetChromStart() int
	GetChromEnd() int
	WriteToFileHandle(io.Writer)
	UpdateCoord(string, int, int) interface{}
}

// GoRead reads Lift interfaces from an input file to a slice. Uses Chans under the hood.
func GoRead(inputFile string) []Lift {
	var answer []Lift
	ch := GoReadToChan(inputFile)
	for i := range ch {
		answer = append(answer, i)
	}
	return answer
}

// GoReadToChan reads Lift interfaces to a channel from an input file.
func GoReadToChan(inputFile string) <-chan Lift {
	answer := make(chan Lift, 1000)
	go ReadToChan(inputFile, answer)
	return answer
}

// ReadToChan reads from a file to send Lift interfaces to a chan<- Lift.
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

// LiftCoordinatesWithChain returns the update coordinates for an input Lift interface based on an input chain.
func LiftCoordinatesWithChain(c chain.Chain, i Lift) (string, int, int) {
	var newStart, newEnd int
	newStart, _ = chain.TPosToQPos(c, i.GetChromStart())
	/* The minus one/plus one handles when a region ends
	with a structural variant and ensures correct placement of the end in the new assembly. */
	newEnd, _ = chain.TPosToQPos(c, i.GetChromEnd()-1)
	newEnd++        //correction for the GetChromEnd -1 on the line above
	if !c.QStrand { //valid bed formats must have start < end. So this correction is made for intervals lifted to the negative strand.
		newStart, newEnd = newEnd, newStart
		newStart += 1 //these lines are corrections for the open/closed interval start and end.
		newEnd += 1
	}
	return c.QName, newStart, newEnd
}

// IntervalSliceToLift casts a slice of intervals to a slice of Lift interfaces.
func IntervalSliceToLift(b []interval.Interval) []Lift {
	var answer []Lift = make([]Lift, len(b))
	for i := range b {
		answer[i] = b[i].(Lift)
	}
	return answer
}

// MatchOverlapLen returns the number of bases shared between two start and end points. Used in MatchProportion.
func MatchOverlapLen(start1 int, end1 int, start2 int, end2 int) int {
	return numbers.Max(0, numbers.Min(end1, end2)-numbers.Max(start1, start2))
}

// MatchProportion returns the proportion of bases in the target and query that can be lifted for a particular interval as a pair of floats (propT, propQ).
func MatchProportion(c chain.Chain, i interval.Interval) (float64, float64) {
	var match, dT, dQ int = 0, 0, 0
	var currPos int = c.TStart //starting with strand +/+ case for now.
	if !c.TStrand {
		log.Fatalf("Please format chain files for lift with the target in the positive strand.")
	}
	for j := 0; j < len(c.Alignment); j++ {
		match += MatchOverlapLen(currPos, currPos+c.Alignment[j].Size, i.GetChromStart(), i.GetChromEnd())
		currPos += c.Alignment[j].Size
		dT += MatchOverlapLen(currPos, currPos+c.Alignment[j].TBases, i.GetChromStart(), i.GetChromEnd())
		if MatchOverlapLen(currPos, currPos+c.Alignment[j].TBases, i.GetChromStart(), i.GetChromEnd()) > 0 {
			/* This handles a special case where both TBases and QBases are non-zero,
			as in competing non-aligned bases are present at the same location in target and query. */
			dQ += c.Alignment[j].QBases
		}
		currPos += c.Alignment[j].TBases
	}
	if match == 0 {
		return 0, 0
	}
	return float64(match) / float64(match+dT), float64(match) / float64(match+dQ)
}

// StrictBorderCheck returns true if the TPos of both the ChromStart and ChromEnd of an interval fall within the chain Size, not TBases.
func StrictBorderCheck(c chain.Chain, i interval.Interval) bool {
	var border bool
	_, border = chain.TPosToQPos(c, i.GetChromStart())
	if !border { //only return if false, otherwise we have to check chromEnd.
		return false
	}
	_, border = chain.TPosToQPos(c, i.GetChromEnd()-1) //interval ranges are open right so we want ChromEnd - 1
	return border
}

//everything below is for lifting with AXT

// refCoordToRefIdx is a helper function for LiftCoordinatesWithAxt which identifies the start and end indices in the Axt reference multifa record that corresponds to the interval to be lifted
func refCoordToRefIdx(a axt.Axt, region interval.Interval) (start, end int) {
	var stopLoop int
	var startFound, endFound bool = false, false
	for i := range a.RSeq {
		if stopLoop >= (region.GetChromStart() - (a.RStart - 1)) {
			startFound = true
			break
		}
		if a.RSeq[i] != dna.Gap {
			stopLoop++
		}
		start++
	}
	end = start
	stopLoop = 0
	for i := start; i < len(a.RSeq); i++ {
		if a.RSeq[i] != dna.Gap {
			stopLoop++
		}
		end++
		if stopLoop >= (region.GetChromEnd() - region.GetChromStart()) {
			endFound = true
			break
		}
	}
	if !startFound || !endFound {
		log.Fatalf("Error in refCoordToRefIdx -- the region to lift could not be found within the axt record.\n")
	}
	return start, end
}

// translateCoord is a helper function for LiftCoordinatesWithAxt which takes the indices of the multifa for the region to be lifted and put them in query coordinates
func translateCoord(a axt.Axt, start, end int) (newStart, newEnd int) {
	newStart = (a.QStart - 1) + dna.CountBasesNoGaps(a.QSeq[0:start])
	newEnd = newStart + dna.CountBasesNoGaps(a.QSeq[start:end])
	return newStart, newEnd
}

// LiftCoordinatesWithAxt takes an Axt record, a lift-compatible interval and chromosome size for the query (in case of a minus strand alignment), and lifts the interval from the reference
// coordinates to the query coordinates. The interval must be completely contained within the reference Axt region. The output is chromosome name, start (0-based) and end coordinates
// of the lifted interval
func LiftCoordinatesWithAxt(a axt.Axt, region Lift, Qsize int) (chrom string, start, end int) {
	if !checkCompatability(a, region) {
		log.Fatalf("The interval you are trying to lift is not entirely within the axt refernce coordinates.")
	}
	chrom = a.QName
	refStart, refEnd := refCoordToRefIdx(a, region)
	newStart, newEnd := translateCoord(a, refStart, refEnd)
	start, end = newStart, newEnd
	if !a.QStrandPos {
		tmpStart := Qsize - end
		tmpEnd := Qsize - start + 1
		start, end = tmpStart, tmpEnd
	}
	return chrom, start, end
}

// checkCompatability is a helper function for LiftCoordinatesWithAxt which determines if the interval is safe to lift over given the Axt alignment
func checkCompatability(a axt.Axt, region interval.Interval) bool {
	if a.RName == region.GetChrom() && a.RStart <= (region.GetChromStart()+1) && a.REnd >= region.GetChromEnd() {
		//The region falls completely within the axt region. We can safely lift.
		return true
	}
	return false
}

func AxtPercentIdentityInInterval(a axt.Axt, refInterval interval.Interval) float64 {
	if !checkCompatability(a, refInterval) {
		log.Fatalf("The interval you are trying to assay is not entirely within the axt refernce coordinates.")
	}
	idxStart, idxEnd := refCoordToRefIdx(a, refInterval)
	return percentIdentity(a, idxStart, idxEnd)
}

func percentIdentity(a axt.Axt, idxStart, idxEnd int) float64 {
	var c int
	for i := idxStart; i < idxEnd; i++ {
		if a.RSeq[i] == a.QSeq[i] {
			c++
		}
	}
	return (float64(c) / float64(idxEnd-idxStart)) * 100
}
