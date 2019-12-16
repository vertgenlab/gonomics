package cigar

import ()

func AddCigar(cigs []*Cigar, newCig *Cigar) []*Cigar {
	if len(cigs) == 0 {
		cigs = append(cigs, newCig)
	} else if cigs[len(cigs)-1].Op == newCig.Op {
		cigs[len(cigs)-1].RunLength += newCig.RunLength
	} else {
		cigs = append(cigs, newCig)
	}
	return cigs
}

func CatCigar(cigs []*Cigar, newCigs []*Cigar) []*Cigar {
	if len(newCigs) == 0 {
		return cigs
	} else if len(cigs) == 0 {
		return newCigs
	} else {
		cigs = AddCigar(cigs, newCigs[0])
		cigs = append(cigs, newCigs[1:]...)
		return cigs
	}
}
