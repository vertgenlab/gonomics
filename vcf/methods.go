package vcf

import (
	"github.com/vertgenlab/gonomics/dna"
	"github.com/vertgenlab/gonomics/fileio"
)

// Current methods satisfy requirements for the following interfaces:
// bed.BedLike

func (v *Vcf) GetChrom() string {
	return v.Chr
}

func (v *Vcf) SetExclude() {
	v.Chr = "EXCLUDE"
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
func (v *Vcf) GetChromStart() int {
	refBases := dna.StringToBases(v.Ref)
	if len(refBases) == 1 {
		return int(v.Pos - 1)
	} else {
		return int(v.Pos)
	}
}

func (v *Vcf) GetChromEnd() int {
	refBases := dna.StringToBases(v.Ref)
	if len(refBases) == 1 {
		return int(v.Pos)
	} else {
		return int(v.Pos) + len(refBases) - 1
	}
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
	Write(file, v)
}

func (v *Vcf) WriteToFileHandle(file *fileio.EasyWriter) {
	WriteVcf(file, v)
}

func (v *Vcf) NextRealRecord(file *fileio.EasyReader) bool {
	var done bool
	var next *Vcf
	for next == nil && !done {
		next, done = NextVcf(file)
	}
	if done {
		return true
	}
	*v = *next
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
