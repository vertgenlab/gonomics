package sam

import "github.com/vertgenlab/gonomics/cigar"

func FilterByQuality(a Sam, filter int) bool {
	return int(a.MapQ) >= filter
}

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
