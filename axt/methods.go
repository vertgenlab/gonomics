package axt

// Current methods satisfy requirements for the following interfaces:
// bed.BedLike

func (a *Axt) GetChr() string {
	return a.RName
}

// to conform with bed standards the startpos will be zero base and the endpos will be 1 base

func (a *Axt) GetStart() int {
	return int(a.RStart - 1)
}

func (a *Axt) GetEnd() int {
	return int(a.REnd)
}
