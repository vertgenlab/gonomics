package cigar

import (
//"github.com/vertgenlab/gonomics/dna"
//"github.com/vertgenlab/gonomics/sam"
//"github.com/vertgenlab/gonomics/fasta"
)

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

func CountCigar(cig *Cigar) (int64, int64) {
	var refIdx, altIdx int64 = 0, 0
	if ConsumesReference(cig.Op) {
		refIdx += cig.RunLength
	}
	if ConsumesQuery(cig.Op) {
		altIdx += cig.RunLength
	}
	return refIdx, altIdx
}

func CountIndexCigar(cig []*Cigar) (int64, int64) {
	var refIdx, altIdx int64 = 0, 0
	for i := 0; i < len(cig);i++ {
		if ConsumesReference(cig[i].Op) {
			refIdx += cig[i].RunLength
		}
		if ConsumesQuery(cig[i].Op) {
			altIdx += cig[i].RunLength
		}
	}
	return refIdx, altIdx
}

func UpdateIndices(ref int64, alt int64, cig *Cigar) {
	newRef, newAlt := CountCigar(cig)
	ref+= newRef
	alt+= newAlt
}

func UpdateIndexSlice(ref int64, alt int64, cigs []*Cigar) {
	for _, cig := range cigs {
		UpdateIndices(ref, alt, cig)
	}
}

//TODO: finish riley's voting matrix to handle indels better
/*
type Allele struct {
	//zero based
	RefPos int64
	RefSeq []dna.Base
	AltSeq []dna.Base
}

func SamToAllele(samRecord *sam.SamAln, genome map[string]*[]dna.Base) []*Allele {
	var answer []*Allele
	var refIdx = samRecord.Pos - 1
	var queryIdx = 0
	var j int64
	var chrom *[]dna.Base = genome[samRecord.RName]

	for i := 0; i < len(samRecord.Cigar); i++ {
		switch samRecord.Cigar[i].Op {

		case 'M':

			for j = 0; j < samRecord.Cigar[i].RunLength; j++ {

				curr := Allele{RefPos: refIdx+j, RefSeq: (*chrom)[refIdx+j:refIdx]}
				answer = append(answer, &curr)
			}

		}

	}

}
/*
func ParseCigarFromSam(samFilename string, fileName string, v []*vcf.Vcf) {
	var aln *sam.SamAln = nil
	var done bool = false
	//var err error
	//sam file to read
	//var numDiscardReads int64
	//var code uint64
	var loc Location
	for aln, done = sam.NextAlignment(samFile); done != true; aln, done = sam.NextAlignment(samFile) {

		fSnpCount, mSnpCount = 0, 0
		//use 1-based
		currRefPos = aln.Pos-1
		currQueryPos = 0

		//faIdx = chroms[aln.RName].Order

		for j = 0; j < len(aln.Cigar); j++ {
			switch aln.Cigar[j].Op {
			case 'S':
				//currRefPos += aln.Cigar[j].RunLength
				currQueryPos += aln.Cigar[j].RunLength
			case 'I':
				currQueryPos += aln.Cigar[j].RunLength
			case 'D':
				currRefPos += aln.Cigar[j].RunLength
			case 'M':
				for k = 0; k < aln.Cigar[j].RunLength; k++ {
					//code = chromAndPosToNumber(faIdx, currRefPos+k)
					loc = Location{Chr: aln.RName, Pos: currRefPos+k}
					_, ok = aMap[loc]
					if ok {
						//fmt.Println("Ref: ", dna.BasesToString(fMap[faIdx|currRefPos+k]), "Alt: ", dna.BasesToString(mMap[faIdx|currRefPos+k]))
						if aln.Seq[currQueryPos+k] == aMap[loc].Ref {
							fSnpCount++
						}
						if aln.Seq[currQueryPos+k] == aMap[loc].Alt {
							mSnpCount++
						}
					}
				}
				currRefPos += aln.Cigar[j].RunLength
				currQueryPos += aln.Cigar[j].RunLength
			}
		}
		if fSnpCount > mSnpCount {
			//fresh = append(fresh, aln)
			sam.WriteAlnToFileHandle(fresh, aln)
		}
		if mSnpCount > fSnpCount {
			//marine = append(marine, aln)
			sam.WriteAlnToFileHandle(marine, aln)
		}
	}
}*/
