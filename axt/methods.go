package axt

// Current methods satisfy requirements for the following interfaces:
// bed.BedLike

func (a *Axt) GetChrom() string {
	return a.RName
}

// to conform with bed standards the startpos will be zero base and the endpos will be 1 base

func (a *Axt) GetChromStart() int {
	return int(a.RStart - 1)
}

func (a *Axt) GetChromEnd() int {
	return int(a.REnd)
}
