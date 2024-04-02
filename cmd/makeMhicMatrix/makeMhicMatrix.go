package main

import (
	"flag"
	"fmt"
	"github.com/vertgenlab/gonomics/chromInfo"
	"github.com/vertgenlab/gonomics/exception"
	"github.com/vertgenlab/gonomics/fileio"
	"github.com/vertgenlab/gonomics/mHiC"
	"github.com/vertgenlab/gonomics/matrix"
	"github.com/vertgenlab/gonomics/numbers"
	"log"
	"math"
	"sort"
)

type settings struct {
	inFile     string
	outDir     string
	chromSizes string
	resolution int
	normMethod string
}

func setToOne(m matrix.Matrix) {
	for i := range m {
		for j := range m[i] {
			m[i][j] = 1
		}
	}
}

func valueSubtractMatrix(m matrix.Matrix, val float64) matrix.Matrix {
	n := matrix.Initilize(matrix.Shape(m))
	for i := range m {
		for j := range m[i] {
			n[i][j] = val - m[i][j]
		}
	}
	return n
}

func scaleMatrixBySum(m matrix.Matrix, sum float64) matrix.Matrix {
	n, _ := matrix.Shape(m)
	scalar := sum / (2.0 * float64(n))
	return matrix.Scale(m, scalar, '*')
}

func krNorm(rawMat matrix.Matrix) matrix.Matrix {
	var k int
	var innertol float64
	var Z, p matrix.Matrix
	nRows, _ := matrix.Shape(rawMat)
	onesMatrix := matrix.Initilize(nRows, 1)
	setToOne(onesMatrix)
	Delta := 3.0
	delta := 0.1
	g := 0.9
	etamax := 0.1
	eta := 0.1
	tol := 1e-6
	stopTol := tol * 0.5
	vector := onesMatrix
	rt := math.Pow(tol, 2) //1e-12
	v := matrix.DotProduct(rawMat, vector)
	rk := valueSubtractMatrix(v, 1)
	rhoKm1 := matrix.DotProduct(matrix.Transpose(rk), rk)[0][0]
	rhoKm2 := rhoKm1
	rout := rhoKm1
	rold := rhoKm1
	MVP := 0
	i := 0

	for rout > rt {
		i++
		if i > 30 {
			break
		}
		k = 0
		y := onesMatrix
		innertol = numbers.Max(math.Pow(eta, 2)*rout, rt) // rt = 1e-12
		for rhoKm1 > innertol {
			k++
			if k == 1 {
				//this happens first
				Z = matrix.Combine(rk, v, '/')
				p = Z
				rhoKm1 = matrix.DotProduct(matrix.Transpose(rk), Z)[0][0]
			} else {
				beta := rhoKm1 / rhoKm2
				p = matrix.Combine(Z, matrix.Scale(p, beta, '*'), '+')
			}
			if k > 10 {
				break
			}
			w := matrix.Combine(matrix.Combine(vector, matrix.DotProduct(rawMat, matrix.Combine(vector, p, '*')), '*'), matrix.Combine(v, p, '*'), '+')
			alpha := rhoKm1 / matrix.DotProduct(matrix.Transpose(p), w)[0][0]
			ap := matrix.Scale(p, alpha, '*')
			ynew := matrix.Combine(y, ap, '+')
			if matrix.Minimum(ynew) <= delta {
				ind := matrix.Where(ap, "<", 0.0)[0]
				gamma := matrix.Minimum(matrix.Combine(valueSubtractMatrix(matrix.SubsetRows(y, ind), delta), matrix.SubsetRows(ap, ind), '/'))
				y = matrix.Combine(y, matrix.Scale(ap, gamma, '*'), '+')
				break
			}
			if matrix.Maximum(ynew) >= Delta {
				ind := matrix.Where(ynew, ">", Delta)[0]
				gamma := matrix.Minimum(matrix.Combine(valueSubtractMatrix(matrix.SubsetRows(y, ind), Delta), matrix.SubsetRows(ap, ind), '/'))
				y = matrix.Combine(y, matrix.Scale(ap, gamma, '*'), '+')
				break
			}
			y = ynew
			rk = matrix.Combine(rk, matrix.Scale(w, alpha, '*'), '-')
			rhoKm2 = rhoKm1
			Z = matrix.Combine(rk, v, '/')
			rhoKm1 = matrix.DotProduct(matrix.Transpose(rk), rk)[0][0]
		}
		vector = matrix.Combine(vector, y, '*')
		v = matrix.Combine(vector, matrix.DotProduct(rawMat, vector), '*')
		rk = valueSubtractMatrix(v, 1.0)
		rhoKm1 = matrix.DotProduct(matrix.Transpose(rk), rk)[0][0]
		rout = rhoKm1
		MVP += k + 1

		rat := rout / rold
		rold = rout
		resNorm := math.Pow(rat, 0.5)
		eta_o := eta
		eta = g * rat
		if math.Pow(g*eta_o, 2) > 0.1 {
			eta = numbers.Max(eta, math.Pow(g*eta_o, 2.0))
		}
		eta = numbers.Max(numbers.Min(eta, etamax), stopTol/resNorm)
	}
	return vector
}

func oneDividedByMatrix(m matrix.Matrix) matrix.Matrix {
	var j int
	o := matrix.Initilize(matrix.Shape(m))
	for i := range m {
		for j = range m[i] {
			if m[i][j] == 0 {
				o[i][j] = 0
			} else {
				o[i][j] = 1.0 / m[i][j]
			}
		}
	}
	return o
}

func enumerate(s []float64) matrix.Matrix {
	m := matrix.Initilize(len(s), 2)
	for i := range s {
		m[i] = []float64{float64(i), s[i]}
	}
	return m
}

func sortSums(m matrix.Matrix) {
	sort.Slice(m, func(i, j int) bool {
		return m[i][1] < m[j][1]
	})
}

func sortIdx(m matrix.Matrix) {
	sort.Slice(m, func(i, j int) bool {
		return m[i][0] < m[j][0]
	})
}

func zeroRows(m matrix.Matrix) matrix.Matrix {
	var n matrix.Matrix
	//get the indexes of the rows/columns to remove
	sums := matrix.AllSums(m)
	for i := range sums[0] {
		if sums[0][i] == 0 {
			continue
		}
		n = append(n, m[i])
	}
	n2 := dropCols(n, sums[1])
	return n2
}

func dropCols(m matrix.Matrix, colSums []float64) matrix.Matrix {
	var n matrix.Matrix = make([][]float64, len(m))
	var dropIdx []int
	for i := range colSums {
		if colSums[i] != 0 {
			continue
		}
		dropIdx = append(dropIdx, i)
	}
	if len(dropIdx) == 0 {
		return m
	}
	for j := range m {
		for i := range dropIdx {
			switch {
			case i == 0:
				n[j] = append(n[j], m[j][0:dropIdx[i]]...)
			case i > 0 && i < len(dropIdx):
				n[j] = append(n[j], m[j][dropIdx[i-1]+1:dropIdx[i]]...)
			}
		}
		n[j] = append(n[j], m[j][dropIdx[len(dropIdx)-1]+1:]...)
	}
	return n
}

func computeBiasVector(s settings, m matrix.Matrix, chrom string, biasOut *fileio.EasyWriter) {
	x := oneDividedByMatrix(m)
	avg := matrix.Average(x)
	fmt.Println(chrom, avg)
	biasVec := matrix.Scale(x, avg, '/')
	for i := range biasVec {
		fileio.WriteToFileHandle(biasOut, fmt.Sprintf("%s\t%d\t%f", chrom, i*s.resolution, biasVec[i][0]))
	}
}

func initializeChromMatix(chromSizes map[string]chromInfo.ChromInfo, chrom string, resolution int) matrix.Matrix {
	numBins := ((chromSizes[chrom].Size) / resolution) + 1
	mat := matrix.Initilize(numBins, numBins)
	return mat
}

func writeMatrix(chrom, name string, mat matrix.Matrix, s settings) {
	out := fileio.EasyCreate(fmt.Sprintf("%s/%s.%s.txt", s.outDir, chrom, name))
	matrix.Write(out, mat, "\t")
	err := out.Close()
	exception.PanicOnErr(err)
}

func writeUniContactsFile(chrom string, m matrix.Matrix, o *fileio.EasyWriter, s settings) {
	var j int
	for i := len(m) - 1; i >= 0; i-- {
		for j = len(m[i]) - 1; j >= i; j-- {
			if m[i][j] != 0 {
				fileio.WriteToFileHandle(o, fmt.Sprintf("%s\t%d\t%d\t%d", chrom, i*s.resolution, j*s.resolution, int(m[i][j])))
			}
		}
	}
}

func krNormalizeMatrix(s settings, m matrix.Matrix, chrom string, biasOut *fileio.EasyWriter) {
	n := zeroRows(m)
	result := krNorm(n)
	computeBiasVector(s, result, chrom, biasOut)
	x := matrix.Diags(matrix.Flatten(result)[0])
	normMat := matrix.DotProduct(x, matrix.DotProduct(n, x))
	finalMatrix := scaleMatrixBySum(normMat, matrix.Sum(n))

	o := fileio.EasyCreate(fmt.Sprintf("%s/%s.krNormMulti.txt", s.outDir, chrom))
	matrix.Write(o, finalMatrix, "\t")
	err := o.Close()
	exception.PanicOnErr(err)
}

func makeMhicMatrix(s settings) {
	var currChrom string
	var matMulti, matUni matrix.Matrix
	uniContactsOut := fileio.EasyCreate(fmt.Sprintf("%s/uniContacts.txt", s.outDir))
	biasOut := fileio.EasyCreate(fmt.Sprintf("%s/biasFile.txt", s.outDir))
	chromSizes := chromInfo.ReadToMap(s.chromSizes)
	interactionChan := mHiC.GoReadInteractionToChan(s.inFile)
	for i := range interactionChan {
		if i.InterChrom {
			continue
		}
		if i.Chrom1 != currChrom {
			if currChrom != "" {
				writeMatrix(currChrom, "multiRaw", matMulti, s)
				krNormalizeMatrix(s, matUni, currChrom, biasOut)
				writeMatrix(currChrom, "uniRaw", matUni, s)
				writeUniContactsFile(currChrom, matUni, uniContactsOut, s)
			}
			currChrom = i.Chrom1
			matMulti = initializeChromMatix(chromSizes, currChrom, s.resolution)
			matUni = initializeChromMatix(chromSizes, currChrom, s.resolution)
		}
		if i.Uni {
			matUni[i.Bin1][i.Bin2]++
			matUni[i.Bin2][i.Bin1]++
		}
		matMulti[i.Bin1][i.Bin2]++
		matMulti[i.Bin2][i.Bin1]++
	}
	err := uniContactsOut.Close()
	exception.PanicOnErr(err)
	err = biasOut.Close()
	exception.PanicOnErr(err)
}

func usage() {
	fmt.Print("makeMhicMatrix -- create interaction matricies from a validPair file.\n" +
		"Input validPair file must be sorted by chrom" +
		"Usage:\n" +
		"makeMhicMatrix [options] validPairFile chromSizesFile\n" +
		"options:\n")
	flag.PrintDefaults()
}

func main() {
	var outDir *string = flag.String("outDir", ".", "Specify a directory to write matrix files to. Default is current directory.")
	var resolution *int = flag.Int("resolution", 10000, "Specify a resolution for binning. Must be the same resolution as used in characterizePairs.go. Default: 10000")
	var normMethod *string = flag.String("normMethod", "", "Specify a matrix normalization method. Options: KR, ICE. Default: NONE")

	flag.Usage = usage
	flag.Parse()
	var expectedNumArgs int = 2
	if len(flag.Args()) != expectedNumArgs {
		flag.Usage()
		log.Fatalf("Expecting %d args, got %d", expectedNumArgs, len(flag.Args()))
	}

	s := settings{
		resolution: *resolution,
		outDir:     *outDir,
		inFile:     flag.Arg(0),
		chromSizes: flag.Arg(1),
		normMethod: *normMethod,
	}
	makeMhicMatrix(s)
}
