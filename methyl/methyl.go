package main

import (
	"github.com/vertgenlab/gonomics/cigar"
	"github.com/vertgenlab/gonomics/dna"
	"github.com/vertgenlab/gonomics/fasta"
	"github.com/vertgenlab/gonomics/fileio"
	"github.com/vertgenlab/gonomics/sam"
)

type MethylCount struct {
	ReadID 	string
	Chr 	string
	Pos 	int
	CG 	int
	RefCG 	int
	CHH 	int
	RefCHH 	int
	CHG 	int
	RefCHG	int
}

func CountCysteine(samFilename string, refFilename string) []*MethylCount {
	answer := make([]*MethylCount, 0)
	samFile := fileio.EasyOpen(samFilename)
	defer samFile.Close()
	var done = false
	var aln *sam.SamAln
	var current *MethylCount

	refFasta := fasta.Read(refFilename)
	fasta.AllToUpper(refFasta)

	ref := fasta.FastaToMap(refFasta)

	for aln, done = sam.NextAlignment(samFile); done != true; aln, done = sam.NextAlignment(samFile) {
		if aln.Cigar[0].Op == '*' {
			continue
		}

		var seqIdx int64 = 0
		var refIdx = aln.Pos - 1

		current = &MethylCount{aln.QName, aln.RName, int(aln.Pos), 0, 0, 0, 0, 0, 0}

		for i := 0; i < len(aln.Cigar); i++ {

			if aln.Cigar[i].Op == 'D' {
				refIdx += aln.Cigar[i].RunLength
			} else if aln.Cigar[i].Op == 'I' {
				seqIdx += aln.Cigar[i].RunLength
			} else if cigar.CigarConsumesReference(*aln.Cigar[i]) {

				for j := 0; int64(j) < aln.Cigar[i].RunLength; j++ {

					if int(seqIdx) > len(aln.Seq) - 3 {
						break
					}


					if ref[aln.RName][refIdx] == dna.C {

						switch {
						case ref[aln.RName][refIdx+1] == dna.G:
							current.RefCG++

						case ref[aln.RName][refIdx+1] != dna.G && ref[aln.RName][refIdx+2] == dna.G:
							current.RefCHG++

						case ref[aln.RName][refIdx+1] != dna.G && ref[aln.RName][refIdx+2] != dna.G:
							current.RefCHH++
						}

						if aln.Seq[seqIdx] == dna.C {
							switch {
							case aln.Seq[seqIdx+1] == dna.G:
								current.CG++

							case aln.Seq[seqIdx+1] != dna.G && aln.Seq[seqIdx+2] == dna.G:
								current.CHG++

							case aln.Seq[seqIdx+1] != dna.G && aln.Seq[seqIdx+2] != dna.G:
								current.CHH++
							}
						}
					}
					seqIdx++
					refIdx++
				}
			} else if aln.Cigar[i].Op != 'H' {
				seqIdx += aln.Cigar[i].RunLength
			}
		}
		answer = append(answer, current)
	}
	return answer
}
