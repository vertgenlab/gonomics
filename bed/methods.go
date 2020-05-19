package bed

// Current methods satisfy requirements for the following interfaces:
// bed.BedLike

func (b *Bed) GetChr() string {
	return b.Chrom
}

func (b *Bed) GetStart() int {
	return int(b.ChromStart)
}

func (b *Bed) GetEnd() int {
	return int(b.ChromEnd)
}
