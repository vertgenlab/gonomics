package net

// GetChrom returns the target chromosome name. This helps implement the Interval interface
func (n Net) GetChrom() string {
	return n.TName
}

// GetChromStart returns the target chromosome start. This helps implement the Interval interface
func (n Net) GetChromStart() int {
	return n.TStart
}

// GetChromEnd returns the target chromosome start. This helps implement the Interval interface
func (n Net) GetChromEnd() int {
	return n.TStart + n.TSize
}
