package convert

import (
	"log"
)

func EnsemblToUCSCGeneral(in string) string {
	//TODO: can we make in a number
	//then
	return "chr" + in
}

//TODO UCSCToEnsemblGeneral

func EnsemblToUCSC(in string) string {
	switch in {
	case "1":
		return "chr1"
	case "2":
		return "chr2"
	case "3":
		return "chr3"
	case "4":
		return "chr4"
	case "5":
		return "chr5"
	case "6":
		return "chr6"
	case "7":
		return "chr7"
	case "8":
		return "chr8"
	case "9":
		return "chr9"
	case "10":
		return "chr10"
	case "11":
		return "chr11"
	case "12":
		return "chr12"
	case "13":
		return "chr13"
	case "14":
		return "chr14"
	case "15":
		return "chr15"
	case "16":
		return "chr16"
	case "17":
		return "chr17"
	case "18":
		return "chr18"
	case "19":
		return "chr19"
	case "20":
		return "chr20"
	case "21":
		return "chr21"
	case "22":
		return "chr22"
	case "X":
		return "chrX"
	case "Y":
		return "chrY"
	}
	log.Fatalf("chr: %s not found.", in)
	return in
}

func UCSCToEnsembl(in string) string {
	switch in {
	case "chr1":
		return "1"
	case "chr2":
		return "2"
	case "chr3":
		return "3"
	case "chr4":
		return "4"
	case "chr5":
		return "5"
	case "chr6":
		return "6"
	case "chr7":
		return "7"
	case "chr8":
		return "8"
	case "chr9":
		return "9"
	case "chr10":
		return "10"
	case "chr11":
		return "11"
	case "chr12":
		return "12"
	case "chr13":
		return "13"
	case "chr14":
		return "14"
	case "chr15":
		return "15"
	case "chr16":
		return "16"
	case "chr17":
		return "17"
	case "chr18":
		return "18"
	case "chr19":
		return "19"
	case "chr20":
		return "20"
	case "chr21":
		return "21"
	case "chr22":
		return "22"
	case "chrX":
		return "X"
	case "chrY":
		return "Y"
	}
	log.Fatalf("chr: %s not found.", in)
	return in
}
