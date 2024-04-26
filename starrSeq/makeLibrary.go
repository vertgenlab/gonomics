package starrSeq

import (
	"fmt"
	"github.com/vertgenlab/gonomics/dna"
	"github.com/vertgenlab/gonomics/exception"
	"github.com/vertgenlab/gonomics/fasta"
	"github.com/vertgenlab/gonomics/fileio"
	"math"
	"strings"
)

type MakeLibSettings struct {
	UpHA        string
	DownHA      string
	Constructs  string
	CS          string
	UmiLen      int
	MinTemp     float64
	MaxTemp     float64
	MaxOligoLen int
	OutFile     string
	SumFile     string
}

func checkTm(fa fasta.Fasta, mid int, minTemp, maxTemp float64) (bool, float64) {
	Tm := dna.MeltingTemp(fa.Seq[mid-10 : mid+10])
	if Tm >= minTemp && Tm <= maxTemp {
		return true, Tm
	}
	return false, Tm
}

func createOligos(seq fasta.Fasta, up, down int) []fasta.Fasta {
	var oligos []fasta.Fasta
	oligo := fasta.Fasta{Name: seq.Name + "_fwd"}
	oligo.Seq = seq.Seq[:up]
	oligos = append(oligos, oligo)
	oligo = fasta.Fasta{Name: seq.Name + "_rev"}
	oligo.Seq = dna.ReverseComplementAndCopy(seq.Seq[down:])
	oligos = append(oligos, oligo)
	return oligos
}

func extendOverlap(bestMid int, minTemp float64, maxOligoSize int, construct fasta.Fasta) (Tm float64, up int, down int) {
	up, down = bestMid+10, bestMid-10
	var oob bool = false

	for i := 0; i < 5; i++ {
		up++
		if up > maxOligoSize {
			up--
			oob = true
			break
		}
		Tm = dna.MeltingTemp(construct.Seq[down:up])
		if Tm >= minTemp {
			return Tm, up, down
		}
	}
	if oob {
		for i := 0; i < 5; i++ {
			down--
			if down < len(construct.Seq)-maxOligoSize {
				down++
				break
			}
			if up-down >= 25 {
				break
			}
			Tm = dna.MeltingTemp(construct.Seq[down:up])
			if Tm >= minTemp {
				return Tm, up, down
			}
		}
	}
	return dna.MeltingTemp(construct.Seq[down:up]), up, down
}

func optimizeMeltingTemps(constructs []fasta.Fasta, minTemp, maxTemp float64, maxOligoSize int, summaryFile string) []fasta.Fasta {
	var mid, c, bestMid, up, down int
	var pass, odd, inRange bool = false, true, false
	var oligos []fasta.Fasta
	var Tm, bestDeltaTm, newTm float64
	var o *fileio.EasyWriter

	if summaryFile != "" {
		o = fileio.EasyCreate(summaryFile)
		fileio.WriteToFileHandle(o, "construct\tTm\toverlapBases\toligoSizeF\toligoSizeR\twithinRange\toverlapSize\tconstructSize")
	}

	for i := range constructs {
		mid = len(constructs[i].Seq) / 2
		inRange, Tm = checkTm(constructs[i], mid, minTemp, maxTemp)

		if inRange {
			oligos = append(oligos, createOligos(constructs[i], mid+10, mid-10)...)
			if summaryFile != "" {
				fileio.WriteToFileHandle(o, fmt.Sprintf("%s\t%f\t%d-%d\t%d\t%d\ttrue\t20\t%d", constructs[i].Name, Tm, mid-10, mid+10, mid+10, len(constructs[i].Seq)-mid+10, len(constructs[i].Seq)))
			}
			continue
		}

		bestMid = mid
		bestDeltaTm = Tm
		c = 0
		pass = false

		for !pass {
			c++
			if odd {
				mid += c
				odd = false
			} else {
				mid -= c
				odd = true
			}
			if (!odd && mid+10 > maxOligoSize) || (odd && mid-10 < len(constructs[i].Seq)-maxOligoSize) {
				if bestDeltaTm < minTemp {
					newTm, up, down = extendOverlap(bestMid, minTemp, maxOligoSize, constructs[i])
					if newTm > bestDeltaTm {
						if newTm >= minTemp && newTm <= maxTemp {
							inRange = true
						} else {
							inRange = false
						}
						oligos = append(oligos, createOligos(constructs[i], up, down)...)
						if summaryFile != "" {
							fileio.WriteToFileHandle(o, fmt.Sprintf("%s\t%f\t%d-%d\t%d\t%d\t%t\t%d\t%d", constructs[i].Name, newTm, down, up, up, len(constructs[i].Seq)-down, inRange, up-down, len(constructs[i].Seq)))
						}
						break
					}
				}
				fmt.Println("couldn't find an ideal mid. using this instead: ", bestMid, bestDeltaTm, constructs[i].Name)
				oligos = append(oligos, createOligos(constructs[i], bestMid+10, bestMid-10)...)
				if summaryFile != "" {
					fileio.WriteToFileHandle(o, fmt.Sprintf("%s\t%f\t%d-%d\t%d\t%d\tfalse\t20\t%d", constructs[i].Name, bestDeltaTm, bestMid-10, bestMid+10, bestMid+10, len(constructs[i].Seq)-bestMid+10, len(constructs[i].Seq)))
				}
				break
			}
			inRange, Tm = checkTm(constructs[i], mid, minTemp, maxTemp)
			if inRange {
				oligos = append(oligos, createOligos(constructs[i], mid+10, mid-10)...)
				if summaryFile != "" {
					fileio.WriteToFileHandle(o, fmt.Sprintf("%s\t%f\t%d-%d\t%d\t%d\ttrue\t20\t%d", constructs[i].Name, Tm, mid-10, mid+10, mid+10, len(constructs[i].Seq)-mid+10, len(constructs[i].Seq)))
				}
				pass = true
				continue
			} else {
				if math.Abs(((maxTemp+minTemp)/2)-Tm) < math.Abs(((maxTemp+minTemp)/2)-bestDeltaTm) {
					bestDeltaTm = Tm
					bestMid = mid
				}
			}
		}
	}
	if summaryFile != "" {
		err := o.Close()
		exception.PanicOnErr(err)
	}
	return oligos
}

func MakeLibrary(s MakeLibSettings) {
	var cs, constructs, oligos []fasta.Fasta
	var construct fasta.Fasta
	var umi []dna.Base
	var tmpUMI []string
	var mid int

	testFa := fasta.Read(s.Constructs)
	up := fasta.Read(s.UpHA)
	down := fasta.Read(s.DownHA)

	if s.CS != "" {
		cs = fasta.Read(s.CS)
	}

	if s.UmiLen != 0 {
		for i := 0; i < s.UmiLen; i++ {
			tmpUMI = append(tmpUMI, "N")
		}
		umi = dna.StringToBases(strings.Join(tmpUMI, ""))
	}

	for i := range testFa {
		construct = fasta.Fasta{Name: testFa[i].Name, Seq: up[0].Seq}
		construct.Seq = append(construct.Seq, testFa[i].Seq...)
		construct.Seq = append(construct.Seq, umi...)
		construct.Seq = append(construct.Seq, cs[0].Seq...)
		construct.Seq = append(construct.Seq, down[0].Seq...)
		constructs = append(constructs, construct)
	}
	if s.MinTemp != -1 || s.MaxTemp != -1 {
		oligos = optimizeMeltingTemps(constructs, s.MinTemp, s.MaxTemp, s.MaxOligoLen, s.SumFile)
	} else {
		for i := range constructs {
			mid = len(constructs[i].Seq) / 2
			oligos = append(oligos, createOligos(constructs[i], mid+10, mid-10)...)
		}
	}
	fasta.Write(s.OutFile, oligos)
}
