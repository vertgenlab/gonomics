package obo

// AllAreEqual returns true if two input slices of Obo structs are identical.
func AllAreEqual(a []Obo, b []Obo) bool {
	if len(a) != len(b) {
		return false
	}
	for i := range a {
		if !Equal(a[i], b[i]) {
			return false
		}
	}
	return true
}

// Equal returns true if two input Obo entries have the same struct contents.
func Equal(a Obo, b Obo) bool {
	if a.Id != b.Id {
		return false
	}
	if a.Name != b.Name {
		return false
	}
	if a.NameSpace != b.NameSpace {
		return false
	}
	if a.Def != b.Def {
		return false
	}
	if !isAEqual(a.IsA, b.IsA) {
		return false
	}
	if !stringSliceEqual(a.Synonyms, b.Synonyms) {
		return false
	}
	if !stringSliceEqual(a.XRefs, b.XRefs) {
		return false
	}
	if !stringSliceEqual(a.AltIds, b.AltIds) {
		return false
	}
	if !stringSliceEqual(a.Relationships, b.Relationships) {
		return false
	}
	if !stringSliceEqual(a.Comments, b.Comments) {
		return false
	}
	if !otherFieldsEqual(a.OtherFields, b.OtherFields) {
		return false
	}
	return true
}

func isAEqual(a []IsADescription, b []IsADescription) bool {
	var i, j int
	if len(a) != len(b) {
		return false
	}
	for i = range a {
		if a[i].ParentId != b[i].ParentId {
			return false
		}
		if len(a[i].ParentInfo) != len(b[i].ParentInfo) {
			return false
		}
		for j = range a[i].ParentInfo {
			if a[i].ParentInfo[j] != b[i].ParentInfo[j] {
				return false
			}
		}
	}
	return true
}

func stringSliceEqual(a []string, b []string) bool {
	if len(a) != len(b) {
		return false
	}
	for i := range a {
		if a[i] != b[i] {
			return false
		}
	}
	return true
}

func otherFieldsEqual(a map[string][]string, b map[string][]string) bool {
	if len(a) != len(b) {
		return false
	}
	var foundInMap bool
	for currKey := range a {
		if _, foundInMap = b[currKey]; !foundInMap {
			return false
		}
		if !stringSliceEqual(a[currKey], b[currKey]) {
			return false
		}
	}
	return true
}

// EqualHeader returns true if two input Header structs have the same data.
func EqualHeader(a Header, b Header) bool {
	return stringSliceEqual(a.Text, b.Text)
}
