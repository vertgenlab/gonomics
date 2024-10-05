package gff3

// equal compares two Gff3 structs for equality.
func Equal(a, b Gff3) bool {
	// Compare basic fields
	if a.Id != b.Id || a.Source != b.Source || a.Type != b.Type ||
		a.Start != b.Start || a.End != b.End || a.Score != b.Score ||
		a.Strand != b.Strand || a.Phase != b.Phase {
		return false
	}

	// Compare attributes
	if len(a.Attributes) != len(b.Attributes) {
		return false
	}
	for i := range a.Attributes {
		if a.Attributes[i] != b.Attributes[i] {
			return false
		}
	}

	return true
}
