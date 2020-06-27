package axt

type AxtHeavy struct {
	*Axt
	Chrom      string
	ChromSize  int
	ChromStart int
	ChromEnd   int
}

func (a *AxtHeavy) GetChrom() string {
	return a.Chrom
}
func (a *AxtHeavy) GetChromStart() int {
	return a.ChromStart
}
func (a *AxtHeavy) GetChromEnd() int {
	return a.ChromEnd
}