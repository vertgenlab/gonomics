package axt

type AxtHeavy struct {
	Axt
	ChromSize  int
}

func (a *AxtHeavy) GetChrom() string {
	return a.Axt.GetChrom()
}
func (a *AxtHeavy) GetChromStart() int {
	return a.Axt.GetChromStart()
}
func (a *AxtHeavy) GetChromEnd() int {
	return a.Axt.GetChromEnd()
}
