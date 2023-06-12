package bedpe

// GetChrom returns the chrom of a half of the bedpe struct.
func (b BedPeHalf) GetChrom() string {
	return b.Chrom
}

// GetChromStart returns the starting coordinates of the bed struct.
func (b BedPeHalf) GetChromStart() int {
	return b.ChromStart
}

// GetChromEnd returns the end coordinates of the bed struct.
func (b BedPeHalf) GetChromEnd() int {
	return b.ChromEnd
}
