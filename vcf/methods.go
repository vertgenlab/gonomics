package vcf

import (
	"fmt"
	"github.com/vertgenlab/gonomics/dna"
	"github.com/vertgenlab/gonomics/fileio"
	"io"
	"strings"
)

// String implements the fmt.Stringer interface for easy printing of Vcf with the fmt package.
func (v Vcf) String() string {
	return fmt.Sprintf("%s\t%v\t%s\t%s\t%s\t%v\t%s\t%s\t%s\t%s\n", v.Chr, v.Pos, v.Id, v.Ref, strings.Join(v.Alt, ","), v.Qual, v.Filter, v.Info, v.Format, SamplesToString(v.Samples))
}

// String implements the fmt.Stringer interface for easy printing of Sample with the fmt package.
func (s Sample) String() string {
	return sampleToString(s)
}

func (v Vcf) GetChrom() string {
	return v.Chr
}

// to conform with bed standards the startpos will be zero base and the endpos will be 1 base
// Example:
// bed reg  |----|
// 1-base 1  2  3  4
// 0-base 0  1  2  3
// seq    A  C  G  T
// If we have a vcf of 1 ACG -> A
// The return would be
// START = 1
// END = 3

// vcf is read in as 1-base so subtract 1 from v.pos
// for indels, vcf records the startpos as the base prior to the change
// to find the region actually being changed we need to check if it is indel
// and adjust accordingly
func (v Vcf) GetChromStart() int {
	refBases := dna.StringToBases(v.Ref)
	if len(refBases) == 1 {
		return v.Pos - 1
	} else {
		return v.Pos
	}
}

func (v Vcf) GetChromEnd() int {
	refBases := dna.StringToBases(v.Ref)
	if len(refBases) == 1 {
		return v.Pos
	} else {
		return v.Pos + len(refBases) - 1
	}
}

func (v Vcf) UpdateCoord(c string, start int, end int) interface{} {
	v.Chr = c
	v.Pos = start + 1
	return v
}

type VcfSlice []*Vcf

func (v VcfSlice) Len() int { return len(v) }

func (v VcfSlice) Swap(i, j int) { v[i], v[j] = v[j], v[i] }

func (v *VcfSlice) Push(x interface{}) {
	answer := x.(*Vcf)
	*v = append(*v, answer)
}

func (v *VcfSlice) Pop() interface{} {
	oldQueue := *v
	n := len(oldQueue)
	answer := oldQueue[n-1]
	*v = oldQueue[:n-1]
	return answer
}

func (v VcfSlice) Write(file string) {
	Write(file, convertToNonPtr(v))
}

//TODO remove this function once interval is updated
func convertToNonPtr(v []*Vcf) []Vcf {
	answer := make([]Vcf, len(v))
	for i := range v {
		answer[i] = *v[i]
	}
	return answer
}

func (v Vcf) WriteToFileHandle(file io.Writer) {
	WriteVcf(file, v)
}

func (v *Vcf) NextRealRecord(file *fileio.EasyReader) bool {
	var done bool
	var next Vcf
	for next.Chr == "" && !done {
		next, done = NextVcf(file)
	}
	*v = next
	if done {
		return true
	}
	return done
}

func (v *Vcf) Copy() interface{} {
	var answer *Vcf = new(Vcf)
	*answer = *v
	return answer
}

type ByGenomicCoordinates struct {
	VcfSlice
}

func (g ByGenomicCoordinates) Less(i, j int) bool {
	// First sort criteria is chromosome
	if g.VcfSlice[i].GetChrom() < g.VcfSlice[j].GetChrom() {
		return true
	} else if g.VcfSlice[i].GetChrom() == g.VcfSlice[j].GetChrom() {
		// If chroms are equal then sort by start position
		if g.VcfSlice[i].GetChromStart() < g.VcfSlice[j].GetChromStart() {
			return true
		} else if g.VcfSlice[i].GetChromStart() == g.VcfSlice[j].GetChromStart() {
			// If start positions are equal then the shorter region wins
			if g.VcfSlice[i].GetChromEnd() < g.VcfSlice[j].GetChromEnd() {
				return true
			}
		}
	}
	return false
}
