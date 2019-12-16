package alleles

import (
	"fmt"
	"io/ioutil"
	"stats"
)


type MultiplexAlleleCount struct {
	Sample	string
	Chr		string
	Pos 	int32
	BaseA	int32
	BaseC 	int32
	BaseG 	int32
	BaseT 	int32
	Ins 	int32
	Del 	int32
}

type MultiplexSortedAlleles struct {
	Sample	string
	Chr		string
	Pos 	int32
	BaseA	int32
	BaseC 	int32
	BaseG 	int32
	BaseT 	int32
	Ins 	int32
	Del 	int32
	BkgdA	int32 	//Bkgd element is equal to sum(sample[@].BaseN) - sample[i].BaseN
	BkgdC	int32 	//This will be used as the background when running Fishers Exact test in R
	BkgdG	int32
	BkgdT	int32
	BkgdIns	int32
	BkgdDel	int32
}

/*
type BatchFishersE struct {
	Sample	string
	Chr		string
	Pos 	int32
	BaseA	int32
	BaseC 	int32
	BaseG 	int32
	BaseT 	int32
	Ins 	int32
	Del 	int32
	FetA	*big.Float
	FetC	*big.Float
	FetG	*big.Float
	FetT	*big.Float
	FetIns	*big.Float
	FetDel	*big.Float
}
 */

type BatchZScores struct {
	Sample	string
	Chr 	string
	Pos 	int32
	Coverage	int32
	NormBaseA	float32
	NormBaseC	float32
	NormBaseG	float32
	NormBaseT	float32
	NormIns		float32
	NormDel		float32
	ZscoreA		float32
	ZscoreC		float32
	ZscoreG		float32
	ZscoreT		float32
	ZscoreIns	float32
	ZscoreDel	float32
}

//feed in a directory filled with AlleleCount (.ac) files. Merges the files into a triple nested map (chr[pos][sample] data)
func CreateSampleMap(inDirectory string) map[string]map[int32][]*MultiplexAlleleCount {

	files, _ := ioutil.ReadDir(inDirectory)
	var SampleName string
	var current *MultiplexAlleleCount

	SampleMap := make(map[string]map[int32][]*MultiplexAlleleCount)

	for _, file := range files {

		SampleName = file.Name()
		SamplePath := fmt.Sprintf("%s%s", inDirectory, SampleName)
		AlleleCounts := ReadAlleleCounts(SamplePath)

		var i int
		for i = 0; i < len(AlleleCounts); i++ {

			//if the chromosome has already been added to the matrix, move along
			_, ok := SampleMap[AlleleCounts[i].Chr]

			//if the chromosome is NOT in the matrix, initialize
			if ! ok {
				SampleMap[AlleleCounts[i].Chr] = make(map[int32][]*MultiplexAlleleCount)
			}

			current = &MultiplexAlleleCount{
				Sample:	SampleName,
				Chr:	AlleleCounts[i].Chr,
				Pos:	AlleleCounts[i].Pos,
				BaseA:	AlleleCounts[i].BaseA,
				BaseC:	AlleleCounts[i].BaseC,
				BaseG:	AlleleCounts[i].BaseG,
				BaseT:	AlleleCounts[i].BaseT,
				Ins:	AlleleCounts[i].Ins,
				Del:	AlleleCounts[i].Del}

			SampleMap[AlleleCounts[i].Chr][AlleleCounts[i].Pos] = append(SampleMap[AlleleCounts[i].Chr][AlleleCounts[i].Pos], current)
		}
	}
	return SampleMap
}




func BatchSortAlleles(SampleMap map[string]map[int32][]*MultiplexAlleleCount) map[string]map[int32][]*MultiplexSortedAlleles {

	//TODO: make these parameters configurable
	//var lowerAFBound float32 = 0.005 //0.5%
	//var upperAFBound float32 = 0.2 //20%

	var current *MultiplexSortedAlleles
	var i int

	SortedSampleMap := make(map[string]map[int32][]*MultiplexSortedAlleles)

	for _, pos := range SampleMap {
		for _, data := range pos {

			var BkgdA int32 = 0
			var BkgdC int32 = 0
			var BkgdG int32 = 0
			var BkgdT int32 = 0
			var BkgdIns int32 = 0
			var BkgdDel int32 = 0

			//calculate bkgd for current position
			for i = 0; i < len(data); i++ {
				BkgdA = BkgdA + data[i].BaseA
				BkgdC = BkgdC + data[i].BaseC
				BkgdG = BkgdG + data[i].BaseG
				BkgdT = BkgdT + data[i].BaseT
				BkgdIns = BkgdIns + data[i].Ins
				BkgdDel = BkgdDel + data[i].Del
			}

			for i = 0; i < len(data); i++ {

				//if the chromosome has already been added to the matrix, move along
				_, ok := SortedSampleMap[data[i].Chr]

				//if the chromosome is NOT in the matrix, initialize
				if ! ok {
					SortedSampleMap[data[i].Chr] = make(map[int32][]*MultiplexSortedAlleles)
				}


				current = &MultiplexSortedAlleles{
					Sample:	data[i].Sample,
					Chr: 	data[i].Chr,
					Pos: 	data[i].Pos,
					BaseA:	data[i].BaseA,
					BaseC: 	data[i].BaseC,
					BaseG: 	data[i].BaseG,
					BaseT: 	data[i].BaseT,
					Ins: 	data[i].Ins,
					Del: 	data[i].Del,
					BkgdA:	BkgdA - data[i].BaseA,
					BkgdC:	BkgdC - data[i].BaseC,
					BkgdG:	BkgdG - data[i].BaseG,
					BkgdT:	BkgdT - data[i].BaseT,
					BkgdIns:	BkgdIns - data[i].Ins,
					BkgdDel:	BkgdDel - data[i].Del}

				//TODO: the below line must be removed or put in an if/else statement when configurations are added
				SortedSampleMap[data[i].Chr][data[i].Pos] = append(SortedSampleMap[data[i].Chr][data[i].Pos], current)

				//TODO: add these back in when configurable options are added
				//coverage := data[i].BaseA + data[i].BaseC + data[i].BaseG + data[i].BaseT + data[i].Ins + data[i].Del

				//MafA := float32(data[i].BaseA) / float32(coverage)
				//MafC := float32(data[i].BaseC) / float32(coverage)
				//MafG := float32(data[i].BaseG) / float32(coverage)
				//MafT := float32(data[i].BaseT) / float32(coverage)
				//MafIns := float32(data[i].Ins) / float32(coverage)
				//MafDel := float32(data[i].Del) / float32(coverage)

				//TODO: make every filter below a configurable option
				/*
				if len(data) > 5 { //must be at least 5 samples with coverage to build the distribution
					if coverage > 200 { //must have at least 200 coverage
						switch {
						//putative site must be between AF bounds AND have greater than 5 supporting reads
						case lowerAFBound < MafA && MafA < upperAFBound && data[i].BaseA > 5:
							SortedSampleMap[data[i].Chr][data[i].Pos] = append(SortedSampleMap[data[i].Chr][data[i].Pos], current)

						case lowerAFBound < MafC && MafC < upperAFBound && data[i].BaseC > 5:
							SortedSampleMap[data[i].Chr][data[i].Pos] = append(SortedSampleMap[data[i].Chr][data[i].Pos], current)

						case lowerAFBound < MafG && MafG < upperAFBound && data[i].BaseG > 5:
							SortedSampleMap[data[i].Chr][data[i].Pos] = append(SortedSampleMap[data[i].Chr][data[i].Pos], current)

						case lowerAFBound < MafT && MafT < upperAFBound && data[i].BaseT > 5:
							SortedSampleMap[data[i].Chr][data[i].Pos] = append(SortedSampleMap[data[i].Chr][data[i].Pos], current)

						case lowerAFBound < MafIns && MafIns < upperAFBound && data[i].Ins > 5:
							SortedSampleMap[data[i].Chr][data[i].Pos] = append(SortedSampleMap[data[i].Chr][data[i].Pos], current)

						case lowerAFBound < MafDel && MafDel < upperAFBound && data[i].Del > 5:
							SortedSampleMap[data[i].Chr][data[i].Pos] = append(SortedSampleMap[data[i].Chr][data[i].Pos], current)
						}
					}
				}

				 */
			}
		}
	}
	return SortedSampleMap
}






//Normalizes coverage by calculating the percentage of reads present for each allele, also initializes the BatchAlleleCorrections map so must be done before calculating Z scores
func NormalizeCoverage(SampleMap map[string]map[int32][]*MultiplexAlleleCount) map[string]map[int32][]*BatchZScores {

	var current *BatchZScores
	var coverage float32
	var NormBaseA float32
	var NormBaseC float32
	var NormBaseG float32
	var NormBaseT float32
	var NormIns float32
	var NormDel float32
	var i int

	NormalizedSampleMap := make(map[string]map[int32][]*BatchZScores)

	//go through sample map and normalize coverage
	for _, pos := range SampleMap {
		for _, data := range pos {
			for i = 0; i < len(data); i++ {

				//if the chromosome has already been added to the matrix, move along
				_, ok := NormalizedSampleMap[data[i].Chr]

				//if the chromosome is NOT in the matrix, initialize
				if ! ok {
					NormalizedSampleMap[data[i].Chr] = make(map[int32][]*BatchZScores)
				}

				coverage = float32(data[i].BaseA + data[i].BaseC + data[i].BaseG + data[i].BaseT + data[i].Ins + data[i].Del)

				NormBaseA = float32(data[i].BaseA) / coverage
				NormBaseC = float32(data[i].BaseC) / coverage
				NormBaseG = float32(data[i].BaseG) / coverage
				NormBaseT = float32(data[i].BaseT) / coverage
				NormIns = float32(data[i].Ins) / coverage
				NormDel = float32(data[i].Del) / coverage

				current = &BatchZScores{
					Sample:	data[i].Sample,
					Chr: 	data[i].Chr,
					Pos: 	data[i].Pos,
					Coverage:	int32(coverage),
					NormBaseA:	NormBaseA,
					NormBaseC:	NormBaseC,
					NormBaseG:	NormBaseG,
					NormBaseT:	NormBaseT,
					NormIns:	NormIns,
					NormDel:	NormDel,
					ZscoreA:	0,
					ZscoreC:	0,
					ZscoreG:	0,
					ZscoreT:	0,
					ZscoreIns:	0,
					ZscoreDel:	0}

				NormalizedSampleMap[data[i].Chr][data[i].Pos] = append(NormalizedSampleMap[data[i].Chr][data[i].Pos], current)
			}
		}
	}
	return NormalizedSampleMap
}




//Calculates Z scores for sample map output from NormalizeCoverage function
func CalculateZscores(NormalizedSampleMap map[string]map[int32][]*BatchZScores) map[string]map[int32][]*BatchZScores {
	var ZscoreA float32
	var ZscoreC float32
	var ZscoreG float32
	var ZscoreT float32
	var ZscoreIns float32
	var ZscoreDel float32
	var i int

	//go through sample map and calculate Z scores
	for _, pos := range NormalizedSampleMap {
		for _, data := range pos {

			var Alist []float64 = nil
			var Clist []float64 = nil
			var Glist []float64 = nil
			var Tlist []float64 = nil
			var Inslist []float64 = nil
			var Dellist []float64 = nil

			for i = 0; i < len(data); i++ {
				Alist = append(Alist, float64(data[i].NormBaseA))
				Clist = append(Clist, float64(data[i].NormBaseC))
				Glist = append(Glist, float64(data[i].NormBaseG))
				Tlist = append(Tlist, float64(data[i].NormBaseT))
				Inslist = append(Inslist, float64(data[i].NormIns))
				Dellist = append(Dellist, float64(data[i].NormDel))
			}

			MeanA := stats.Mean(Alist)
			StdDevA := stats.StdDev(Alist)

			MeanC := stats.Mean(Clist)
			StdDevC := stats.StdDev(Clist)

			MeanG := stats.Mean(Glist)
			StdDevG := stats.StdDev(Glist)

			MeanT := stats.Mean(Tlist)
			StdDevT := stats.StdDev(Tlist)

			MeanIns := stats.Mean(Inslist)
			StdDevIns := stats.StdDev(Inslist)

			MeanDel := stats.Mean(Dellist)
			StdDevDel := stats.StdDev(Dellist)


			for i = 0; i < len(data); i++ {

				ZscoreA = (data[i].NormBaseA - float32(MeanA)) / float32(StdDevA)
				ZscoreC = (data[i].NormBaseC - float32(MeanC)) / float32(StdDevC)
				ZscoreG = (data[i].NormBaseG - float32(MeanG)) / float32(StdDevG)
				ZscoreT = (data[i].NormBaseT - float32(MeanT)) / float32(StdDevT)
				ZscoreIns = (data[i].NormIns - float32(MeanIns)) / float32(StdDevIns)
				ZscoreDel = (data[i].NormDel - float32(MeanDel)) / float32(StdDevDel)

				data[i].ZscoreA = ZscoreA
				data[i].ZscoreC = ZscoreC
				data[i].ZscoreG = ZscoreG
				data[i].ZscoreT = ZscoreT
				data[i].ZscoreIns = ZscoreIns
				data[i].ZscoreDel = ZscoreDel

			}
		}
	}
	return NormalizedSampleMap
}


func FilterAF(input map[string]map[int32][]*BatchZScores) map[string]map[int32][]*BatchZScores {
	var i int
	var current *BatchZScores
	var lowerAFBound float32 = 0.01 //1%
	var upperAFBound float32 = 0.20 //20%
	var lowerZBound float32 = 1 //must have z score greater than 1

	FilteredVars := make(map[string]map[int32][]*BatchZScores)

	for _, pos := range input {
		for _, data := range pos {
			for i = 0; i < len(data); i++ {

				//if the chromosome has already been added to the matrix, move along
				_, ok := FilteredVars[data[i].Chr]

				//if the chromosome is NOT in the matrix, initialize
				if ! ok {
					FilteredVars[data[i].Chr] = make(map[int32][]*BatchZScores)
			}

				current = &BatchZScores{
					Sample:	data[i].Sample,
					Chr: 	data[i].Chr,
					Pos: 	data[i].Pos,
					Coverage:	data[i].Coverage,
					NormBaseA:	data[i].NormBaseA,
					NormBaseC:	data[i].NormBaseC,
					NormBaseG:	data[i].NormBaseG,
					NormBaseT:	data[i].NormBaseT,
					NormIns:	data[i].NormIns,
					NormDel:	data[i].NormDel,
					ZscoreA:	data[i].ZscoreA,
					ZscoreC:	data[i].ZscoreC,
					ZscoreG:	data[i].ZscoreG,
					ZscoreT:	data[i].ZscoreT,
					ZscoreIns:	data[i].ZscoreIns,
					ZscoreDel:	data[i].ZscoreDel}

				/*
				_, okk := FilteredVars[data[i].Chr][data[i].Pos]

				if ! okk {
					FilteredVars[data[i].Chr][data[i].Pos] = []*BatchAlleleCorrections{}
				}

				 */

				//add in filtering steps

				if len(data) > 5 { //must be at least 5 samples with coverage to build the distribution
					if data[i].Coverage > 200 { //must have at least 200 coverage
						switch {
						case lowerAFBound < data[i].NormBaseA && data[i].NormBaseA < upperAFBound && data[i].ZscoreA > lowerZBound:
							FilteredVars[data[i].Chr][data[i].Pos] = append(FilteredVars[data[i].Chr][data[i].Pos], current)

						case lowerAFBound < data[i].NormBaseC && data[i].NormBaseC < upperAFBound && data[i].ZscoreC > lowerZBound:
							FilteredVars[data[i].Chr][data[i].Pos] = append(FilteredVars[data[i].Chr][data[i].Pos], current)

						case lowerAFBound < data[i].NormBaseG && data[i].NormBaseG < upperAFBound && data[i].ZscoreG > lowerZBound:
							FilteredVars[data[i].Chr][data[i].Pos] = append(FilteredVars[data[i].Chr][data[i].Pos], current)

						case lowerAFBound < data[i].NormBaseT && data[i].NormBaseT < upperAFBound && data[i].ZscoreT > lowerZBound:
							FilteredVars[data[i].Chr][data[i].Pos] = append(FilteredVars[data[i].Chr][data[i].Pos], current)

						case lowerAFBound < data[i].NormIns && data[i].NormIns < upperAFBound && data[i].ZscoreIns > lowerZBound:
							FilteredVars[data[i].Chr][data[i].Pos] = append(FilteredVars[data[i].Chr][data[i].Pos], current)

						case lowerAFBound < data[i].NormDel && data[i].NormDel < upperAFBound && data[i].ZscoreDel > lowerZBound:
							FilteredVars[data[i].Chr][data[i].Pos] = append(FilteredVars[data[i].Chr][data[i].Pos], current)
						}
					}
				}
			}
		}
	}
	return FilteredVars
}



/*
func CalculateFishersE (input map[string]map[int32][]*MultiplexAlleleCount) map[string]map[int32][]*BatchFishersE {
	var answer = make(map[string]map[int32][]*BatchFishersE)
	var curr *BatchFishersE
	var i int

	//go through sample map and calculate FishersE
	for _, pos := range input {
		for _, data := range pos {

			var SumA int32 = 0
			var SumC int32 = 0
			var SumG int32 = 0
			var SumT int32 = 0
			var SumIns int32 = 0
			var SumDel int32 = 0

			for i = 0; i < len(data); i++ {
				SumA = SumA + data[i].BaseA
				SumC = SumC + data[i].BaseC
				SumG = SumG + data[i].BaseG
				SumT = SumT + data[i].BaseT
				SumIns = SumIns + data[i].Ins
				SumDel = SumDel + data[i].Del
			}

			for i = 0; i < len(data); i++ {

				//if the chromosome has already been added to the matrix, move along
				_, ok := answer[data[i].Chr]

				//if the chromosome is NOT in the matrix, initialize
				if ! ok {
					answer[data[i].Chr] = make(map[int32][]*BatchFishersE)
				}

				ExpMajorAllele := alleles.FindMajorAllele(data[i].BaseA, data[i].BaseC, data[i].BaseG, data[i].BaseT, data[i].Ins, data[i].Del)
				var BkgdMajorAllele int32

				switch ExpMajorAllele {

				case data[i].BaseA:
					BkgdMajorAllele = SumA - data[i].BaseA

				case data[i].BaseC:
					BkgdMajorAllele = SumC - data[i].BaseC

				case data[i].BaseG:
					BkgdMajorAllele = SumG - data[i].BaseG

				case data[i].BaseT:
					BkgdMajorAllele = SumT - data[i].BaseT

				case data[i].Ins:
					BkgdMajorAllele = SumIns - data[i].Ins

				case data[i].Del:
					BkgdMajorAllele = SumDel - data[i].Del
				}

				// row/column order for FET is  r1c1, r1c2, r2c1, r2c2
				FetA := gostats.FishersExactTest(ExpMajorAllele, BkgdMajorAllele, data[i].BaseA, SumA - data[i].BaseA)
				FetC := gostats.FishersExactTest(ExpMajorAllele, BkgdMajorAllele, data[i].BaseC, SumC - data[i].BaseC)
				FetG := gostats.FishersExactTest(ExpMajorAllele, BkgdMajorAllele, data[i].BaseG, SumG - data[i].BaseG)
				FetT := gostats.FishersExactTest(ExpMajorAllele, BkgdMajorAllele, data[i].BaseT, SumT - data[i].BaseT)
				FetIns := gostats.FishersExactTest(ExpMajorAllele, BkgdMajorAllele, data[i].Ins, SumIns - data[i].Ins)
				FetDel := gostats.FishersExactTest(ExpMajorAllele, BkgdMajorAllele, data[i].Del, SumDel - data[i].Del)

				curr = &BatchFishersE{
					Sample:	data[i].Sample,
					Chr: 	data[i].Chr,
					Pos: 	data[i].Pos,
					BaseA:	data[i].BaseA,
					BaseC:	data[i].BaseC,
					BaseG:	data[i].BaseG,
					BaseT:	data[i].BaseT,
					Ins:	data[i].Ins,
					Del:	data[i].Del,
					FetA:	FetA,
					FetC:	FetC,
					FetG:	FetG,
					FetT:	FetT,
					FetIns:	FetIns,
					FetDel:	FetDel}

				answer[data[i].Chr][data[i].Pos] = append(answer[data[i].Chr][data[i].Pos], curr)
			}
		}
	}

	return answer
}

 */





/*
func main() {

	flag.Parse()
	inDirectory := flag.Arg(0)

	var i int
	SampleMap := CreateSampleMap(inDirectory)
	fmt.Printf("Normalizing Coverage\n")
	NormalizedSampleMap := NormalizeCoverage(SampleMap)
	fmt.Printf("Calculating Z Scores\n")
	Zscores := CalculateZscores(NormalizedSampleMap)
	for _, pos := range Zscores {
		for _, data := range pos {
			for i = 0; i < len(data); i++ {
				fmt.Println(data[i])
			}
		}
	}
}
 */