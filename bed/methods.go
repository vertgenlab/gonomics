package bed

// Current methods satisfy requirements for the following interfaces:
// bed.BedLike

func (b *Bed) GetChrom() string {
	return b.Chrom
}

func (b *Bed) GetChromStart() int {
	return int(b.ChromStart)
}

func (b *Bed) GetChromEnd() int {
	return int(b.ChromEnd)
}
