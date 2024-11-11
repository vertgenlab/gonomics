package net

import (
	"fmt"
	"github.com/vertgenlab/gonomics/bed"
	"github.com/vertgenlab/gonomics/fileio"
	"github.com/vertgenlab/gonomics/numbers/parse"
)

// bedAnnotation is struct that helps the ToBed function. It holds data that makes up extended bed fields.
type bedAnnotation struct {
	thickStart  int
	thickEnd    int
	itemRGB     string
	blockCount  int
	blockSizes  []int
	blockStarts []int
}

// ToBed takes in a position-sorted slice of Net and converts it to a slice of Bed. Note: since bed entries
// are written to the slice when their last gap has been processed, the bed slice may not be position sorted. The RGB
// color field corresponds to the level of Net.
func ToBed(nets []Net) []bed.Bed {
	var prevChrom string = nets[0].TName
	var ans, currBed []bed.Bed
	var currAnno []bedAnnotation
	var highestLvl int
	for i := range nets {
		if nets[i].TName != prevChrom {
			ans = append(ans, formatAll(currBed, currAnno)...)
			currAnno = []bedAnnotation{}
			currBed = []bed.Bed{}
			highestLvl = 0
		}
		if nets[i].Class == "fill" {
			if nets[i].Level > highestLvl {
				highestLvl = nets[i].Level
				currBed = append(currBed, bed.Bed{})
				currAnno = append(currAnno, bedAnnotation{})
			} else {
				currAnno[nets[i].Level-1].blockSizes = append(currAnno[nets[i].Level-1].blockSizes, (currBed[nets[i].Level-1].ChromEnd-currBed[nets[i].Level-1].ChromStart)-currAnno[nets[i].Level-1].blockStarts[len(currAnno[nets[i].Level-1].blockStarts)-1])
				currBed[nets[i].Level-1].Annotation = annoToStringSlice(currAnno[nets[i].Level-1], nets[i].Level)
				ans = append(ans, currBed[nets[i].Level-1])
			}
			currBed[nets[i].Level-1] = createBed(nets[i])
			currAnno[nets[i].Level-1] = bedAnnotation{thickStart: nets[i].TStart, thickEnd: nets[i].TStart, blockCount: 1, blockStarts: []int{0}, itemRGB: "0,0,0"}
		} else {
			currAnno[nets[i].Level-1].blockCount++
			currAnno[nets[i].Level-1].blockSizes = append(currAnno[nets[i].Level-1].blockSizes, (nets[i].TStart-currBed[nets[i].Level-1].ChromStart)-currAnno[nets[i].Level-1].blockStarts[len(currAnno[nets[i].Level-1].blockStarts)-1])
			currAnno[nets[i].Level-1].blockStarts = append(currAnno[nets[i].Level-1].blockStarts, (nets[i].TStart-currBed[nets[i].Level-1].ChromStart)+nets[i].TSize)
		}
		prevChrom = nets[i].TName
	}
	ans = append(ans, formatAll(currBed, currAnno)...)
	return ans
}

func createBed(n Net) bed.Bed {
	var b bed.Bed = bed.Bed{
		Chrom:             n.TName,
		ChromStart:        n.TStart,
		ChromEnd:          n.TStart + n.TSize,
		Name:              formatBedName(n),
		Score:             n.Level,
		Strand:            bed.Strand(parse.StrandToRune(n.Orientation)),
		FieldsInitialized: 7,
		Annotation:        []string{},
	}
	return b
}

func formatBedName(n Net) string {
	return fmt.Sprintf("%s_%d_%d", n.QName, n.QStart, n.QStart+n.QSize)
}

func annoToStringSlice(a bedAnnotation, lvl int) []string {
	switch lvl {
	case 1:
		a.itemRGB = "153,204,102"
	case 2:
		a.itemRGB = "246,237,100"
	case 3:
		a.itemRGB = "234,51,35"
	case 4:
		a.itemRGB = "124,199,216"
	case 5:
		a.itemRGB = "182,155,197"
	default:
		a.itemRGB = "0,0,0"
	}
	return []string{fileio.IntToString(a.thickStart), fileio.IntToString(a.thickEnd), a.itemRGB, fileio.IntToString(a.blockCount), fileio.IntSliceToString(a.blockSizes), fileio.IntSliceToString(a.blockStarts)}
}

func formatAll(beds []bed.Bed, anno []bedAnnotation) []bed.Bed {
	for i := range beds {
		if anno[i].blockCount != len(anno[i].blockSizes) {
			anno[i].blockSizes = append(anno[i].blockSizes, (beds[i].ChromEnd-beds[i].ChromStart)-anno[i].blockStarts[len(anno[i].blockStarts)-1])
		}
		beds[i].Annotation = annoToStringSlice(anno[i], i+1)
	}
	return beds
}
