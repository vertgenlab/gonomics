package starrSeq

import (
	"fmt"
	"github.com/vertgenlab/gonomics/dna"
	"github.com/vertgenlab/gonomics/exception"
	"github.com/vertgenlab/gonomics/fasta"
	"github.com/vertgenlab/gonomics/fileio"
	"math"
	"sort"
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
	OutDir      string
	SumFile     string
}

func writeOutPools(pools [][]fasta.Fasta, outDir string) {
	var r int = 65
	for i := 0; i < len(pools); i++ {
		fasta.Write(fmt.Sprintf("%s/oligoLib.pool%c.fa", outDir, r+i), pools[i])
	}
}

// searchConstruct returns true if no similar matches are found in the construct and false if a simiar match is found
func searchConstruct(ref fasta.Fasta, query []dna.Base) bool {
	var testSeq []dna.Base
	for i := 0; i < len(ref.Seq)-len(query); i++ {
		testSeq = ref.Seq[i : i+len(query)]
		if dna.SeqsAreSimilar(testSeq, query, 2) {
			return false
		}
	}
	return true
}

func findSearchSeq(overlap []dna.Base, minTemp float64) []dna.Base {
	for i := 10; i < len(overlap); i++ {
		if dna.MeltingTemp(overlap[0:i]) > minTemp {
			return overlap[0:i]
		}
	}
	return overlap
}

func orderPools(pools [][]fasta.Fasta) [][]int {
	var ans [][]int
	for i := range pools {
		ans = append(ans, []int{i, len(pools[i])})
	}
	sort.Slice(ans, func(i, j int) bool {
		return ans[i][1] < ans[j][1]
	})
	return ans
}

func distributeIntoPool(pools, oligoPools [][]fasta.Fasta, overlap []dna.Base, minTm float64) (int, [][]fasta.Fasta, [][]fasta.Fasta) {
	var match bool
	idx := orderPools(pools)
	searchSeq := findSearchSeq(overlap, minTm)
	for i := range idx {
		match = false
		for j := range pools[idx[i][0]] {
			if searchConstruct(pools[idx[i][0]][j], searchSeq) {
				continue
			} else {
				match = true
				break
			}
		}
		if !match {
			return idx[i][0], pools, oligoPools
		}
	}
	pools = append(pools, []fasta.Fasta{})
	oligoPools = append(oligoPools, []fasta.Fasta{})
	return len(pools) - 1, pools, oligoPools
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

	for i := 0; i < 10; i++ {
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
		for i := 0; i < 10; i++ {
			down--
			if down < len(construct.Seq)-maxOligoSize {
				down++
				break
			}
			if up-down >= 30 {
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

func optimizeMeltingTemps(constructs []fasta.Fasta, minTemp, maxTemp float64, maxOligoSize int, summaryFile, outDir string) [][]fasta.Fasta {
	var mid, c, bestMid, up, down, idx int
	var pass, odd, inRange bool = false, true, false
	var Tm, bestDeltaTm, newTm float64
	var o *fileio.EasyWriter
	var oligoPools, pools [][]fasta.Fasta = [][]fasta.Fasta{}, [][]fasta.Fasta{}

	o = fileio.EasyCreate(outDir + "/" + summaryFile)
	fileio.WriteToFileHandle(o, "construct\tTm\toverlapBases\toligoSizeF\toligoSizeR\twithinRange\toverlapSize\tconstructSize\tGC%\tpool")

	for i := range constructs {
		mid = len(constructs[i].Seq) / 2
		inRange, Tm = checkTm(constructs[i], mid, minTemp, maxTemp)

		if inRange {
			idx, pools, oligoPools = distributeIntoPool(pools, oligoPools, constructs[i].Seq[mid-10:mid+10], minTemp)
			oligoPools[idx] = append(oligoPools[idx], createOligos(constructs[i], mid+10, mid-10)...)
			pools[idx] = append(pools[idx], constructs[i])
			fileio.WriteToFileHandle(o, fmt.Sprintf("%s\t%f\t%d-%d\t%d\t%d\ttrue\t20\t%d\t%f\t%c",
				constructs[i].Name, Tm, mid-10, mid+10, mid+10, len(constructs[i].Seq)-mid+10, len(constructs[i].Seq), dna.GCContent(constructs[i].Seq[mid-10:mid+10]), idx+65))
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
						idx, pools, oligoPools = distributeIntoPool(pools, oligoPools, constructs[i].Seq[down:up], minTemp)
						oligoPools[idx] = append(oligoPools[idx], createOligos(constructs[i], up, down)...)
						pools[idx] = append(pools[idx], constructs[i])
						fileio.WriteToFileHandle(o, fmt.Sprintf("%s\t%f\t%d-%d\t%d\t%d\t%t\t%d\t%d\t%f\t%c",
							constructs[i].Name, newTm, down, up, up, len(constructs[i].Seq)-down, inRange, up-down, len(constructs[i].Seq), dna.GCContent(constructs[i].Seq[down:up]), idx+65))
						break
					}
				}
				idx, pools, oligoPools = distributeIntoPool(pools, oligoPools, constructs[i].Seq[bestMid-10:bestMid+10], minTemp)
				oligoPools[idx] = append(oligoPools[idx], createOligos(constructs[i], bestMid+10, bestMid-10)...)
				pools[idx] = append(pools[idx], constructs[i])
				fileio.WriteToFileHandle(o, fmt.Sprintf("%s\t%f\t%d-%d\t%d\t%d\tfalse\t20\t%d\t%f\t%c",
					constructs[i].Name, bestDeltaTm, bestMid-10, bestMid+10, bestMid+10, len(constructs[i].Seq)-bestMid+10, len(constructs[i].Seq), dna.GCContent(constructs[i].Seq[bestMid-10:bestMid+10]), idx+65))
				break
			}
			inRange, Tm = checkTm(constructs[i], mid, minTemp, maxTemp)
			if inRange {
				idx, pools, oligoPools = distributeIntoPool(pools, oligoPools, constructs[i].Seq[mid-10:mid+10], minTemp)
				oligoPools[idx] = append(oligoPools[idx], createOligos(constructs[i], mid+10, mid-10)...)
				pools[idx] = append(pools[idx], constructs[i])
				fileio.WriteToFileHandle(o, fmt.Sprintf("%s\t%f\t%d-%d\t%d\t%d\ttrue\t20\t%d\t%f\t%c",
					constructs[i].Name, Tm, mid-10, mid+10, mid+10, len(constructs[i].Seq)-mid+10, len(constructs[i].Seq), dna.GCContent(constructs[i].Seq[mid-10:mid+10]), idx+65))
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
	err := o.Close()
	exception.PanicOnErr(err)
	return oligoPools
}

func MakeLibrary(s MakeLibSettings) {
	var cs, constructs, oligos []fasta.Fasta
	var pools [][]fasta.Fasta
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
		pools = optimizeMeltingTemps(constructs, s.MinTemp, s.MaxTemp, s.MaxOligoLen, s.SumFile, s.OutDir)
	} else {
		for i := range constructs {
			mid = len(constructs[i].Seq) / 2
			oligos = append(oligos, createOligos(constructs[i], mid+10, mid-10)...)
		}
		fasta.Write(s.OutDir+"/oligoLib.out.fa", oligos)
	}
	writeOutPools(pools, s.OutDir)
}
