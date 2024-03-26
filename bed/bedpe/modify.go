package bedpe

// ContactsToMidpoints returns a bedpe where the contacts have been collapsed to their midpoints
func ContactsToMidpoints(bps []BedPe) {
	var aMidpoint, bMidpoint int
	for i := range bps {
		aMidpoint = (bps[i].A.ChromStart + bps[i].A.ChromEnd) / 2
		bMidpoint = (bps[i].B.ChromStart + bps[i].B.ChromEnd) / 2

		bps[i].A.ChromStart = aMidpoint
		bps[i].A.ChromEnd = aMidpoint + 1
		bps[i].A.Name = ""
		bps[i].A.Score = 0

		bps[i].B.ChromStart = bMidpoint
		bps[i].B.ChromEnd = bMidpoint + 1
		bps[i].B.Name = ""
		bps[i].B.Score = 0
	}
}
