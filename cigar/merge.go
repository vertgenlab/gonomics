package cigar

import (
	//"fmt"
	
)

func AddCigar(cigs []*Cigar, newCig *Cigar) {
	//fmt.Println(cigs[len(cigs)-1])
	if cigs[len(cigs)-1].Op == newCig.Op {
		cigs[len(cigs)-1].RunLength +=newCig.RunLength
	} else {
		cigs = append(cigs, newCig)
	}
}

func AddAllCigar(cigs []*Cigar, newCigs []*Cigar) {
	for i := 0; i < len(newCigs);i++ {
		if cigs[len(cigs)-1].Op == newCigs[i].Op {
			cigs[len(cigs)-1].RunLength +=newCigs[i].RunLength
		} else {
			cigs = append(cigs, newCigs[i])
		}
	}
}