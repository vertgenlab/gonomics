package cigar

import (
//"fmt"

)

func AddCigar(cigs []*Cigar, newCig *Cigar) []*Cigar {
	//fmt.Println(cigs[len(cigs)-1])
	if cigs[len(cigs)-1].Op == newCig.Op {
		cigs[len(cigs)-1].RunLength += newCig.RunLength
	} else {
		cigs = append(cigs, newCig)
	}
	return cigs
}

func CatCigar(cigs []*Cigar, newCigs []*Cigar) []*Cigar {
	if len(newCigs) == 0 {
		return cigs
	} else if cigs[len(cigs)-1].Op == newCigs[0].Op {
		cigs[len(cigs)-1].RunLength += newCigs[0].RunLength
	} else {
		cigs = append(cigs, newCigs[0])
	}
	cigs = append(cigs, newCigs[1:]...)
	return cigs
}
