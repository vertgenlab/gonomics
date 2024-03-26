package gaf

// RemoveDuplicates filters Gaf records that have the same GoId/DbObjectSymbol pairing.
func RemoveDuplicates(records []Gaf) []Gaf {
	var seenPair = make(map[string]bool) //maps a string of "GeneIdGoTerm" to bool
	var answer []Gaf
	var currPair string
	var foundInMap bool

	for i := range records {
		currPair = records[i].GoId + records[i].DbObjectSymbol
		if _, foundInMap = seenPair[currPair]; !foundInMap {
			seenPair[currPair] = true
			answer = append(answer, records[i])
		}
	}

	return answer
}
