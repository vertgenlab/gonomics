package sam

import "github.com/vertgenlab/gonomics/cigar"

// FilterByQuality takes in a sam record and an integer for MAPQ filtering. Returns true if the read has equal or higher
// MAPQ than the filter. Returns false if the read MAPQ is lower than the filter value
func FilterByQuality(a Sam, filter int) bool {
	return int(a.MapQ) >= filter
}

// FilterByCigar takes in a sam record and a string corresponding to a preset option or cigar string. Returns true if the
// cigar string is an exact match or if it passes the present filter option
func FilterByCigar(a Sam, filter string) bool {
	if filter == "starrSeqIntrons" {
		for _, k := range a.Cigar { //filter out introns longer than 500bp (the length of 1 standard STARR-seq construct)
			if k.Op == 'N' && k.RunLength > 500 {
				return false
			}
		}
	} else {
		if cigar.ToString(a.Cigar) != filter {
			return false
		}
	}
	return true
}
