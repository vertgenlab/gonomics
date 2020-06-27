package axt

type AxtHeavy struct {
	*Axt
	ChromSize  int
}
//TODO: GetChrom, GetStart, GetEnd returns target chr, implement a target and query swap and call same function on target
func (a *AxtHeavy) GetChrom() string {
	return a.RName
}
func (a *AxtHeavy) GetChromStart() int {
	return a.RStart - 1
}
func (a *AxtHeavy) GetChromEnd() int {
	return a.REnd
}
