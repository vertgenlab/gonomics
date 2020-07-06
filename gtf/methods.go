package gtf

func (t *Transcript) GetChrom() string {
	return t.Chr
}

func (t *Transcript) GetChromStart() int {
	t.Start - 1
}