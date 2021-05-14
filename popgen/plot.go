package popgen

import (
	"fmt"
	"github.com/vertgenlab/gonomics/exception"
	"github.com/vertgenlab/gonomics/fileio"
)

//PlotAfsF writes the Allele Frequency F function (AfsSampleDensity) to an output file for downstream visualization.
func PlotAfsF(alpha float64, n int, outFile string, integralError float64) {
	out := fileio.EasyCreate(outFile)
	allN := []int{n}
	binomCache := BuildBinomCache(allN)

	fmt.Fprintf(out, "Frequency\tF\n")

	for i := 1; i < n; i++ {
		fmt.Fprintf(out, "%v\t%e\n", i, AfsSampleDensity(n, i, alpha, binomCache, integralError))
	}
	err := out.Close()
	exception.PanicOnErr(err)
}

//PlotAfsFCareless plots the old version of AfsF, which lacks the "careful" integration. Used for plotting and debugging. TODO: Delete eventually.
func PlotAfsFCareless(alpha float64, n int, outFile string) {
	out := fileio.EasyCreate(outFile)
	fmt.Fprintf(out, "Frequency\tF\n")
	for i := 1; i < n; i++ {
		fmt.Fprintf(out, "%v\t%e\n", i, AfsSampleDensityCareless(n, i, alpha))
	}
	err := out.Close()
	exception.PanicOnErr(err)
}

func PlotAfsPmfAncestral(alpha float64, n int, outFile string, integralError float64) {
	out := fileio.EasyCreate(outFile)

	allN := []int{n}
	binomCache := BuildBinomCache(allN)
	//DEBUG:fmt.Printf("Len BinomCache:%v. Entry at index n is of length: %v.\n", len(binomCache), len(binomCache[n]))

	//a short header is written to the file
	fmt.Fprintf(out, "Frequency\tProbability\n")
	//for each possible allele frequency for a segregating site, we calculate the probability to the file
	for i := 1; i < n; i++ {
		fmt.Fprintf(out, "%v\t%e\n", i, AlleleFrequencyProbabilityAncestralAscertainment(alpha, i, n, 1, binomCache, integralError))
	}

	err := out.Close()
	exception.PanicOnErr(err)
}

func PlotAfsPmfDerived(alpha float64, n int, outFile string, integralError float64) {
	out := fileio.EasyCreate(outFile)

	allN := []int{n}
	binomCache := BuildBinomCache(allN)
	//DEBUG:fmt.Printf("Len BinomCache:%v. Entry at index n is of length: %v.\n", len(binomCache), len(binomCache[n]))

	//a short header is written to the file
	fmt.Fprintf(out, "Frequency\tProbability\n")
	//for each possible allele frequency for a segregating site, we calculate the probability to the file
	for i := 1; i < n; i++ {
		fmt.Fprintf(out, "%v\t%e\n", i, AlleleFrequencyProbabilityDerivedAscertainment(alpha, i, n, 1, binomCache, integralError))
	}

	err := out.Close()
	exception.PanicOnErr(err)
}

//PlotAfsPmf writes the allele frequency probability mass function (AlleleFrequencyProbability) to an output file for downstream visualization.
func PlotAfsPmf(alpha float64, n int, outFile string, integralError float64) {
	out := fileio.EasyCreate(outFile)

	allN := []int{n}
	binomCache := BuildBinomCache(allN)
	//DEBUG:fmt.Printf("Len BinomCache:%v. Entry at index n is of length: %v.\n", len(binomCache), len(binomCache[n]))

	//a short header is written to the file
	fmt.Fprintf(out, "Frequency\tProbability\n")
	//for each possible allele frequency for a segregating site, we calculate the probability to the file
	for i := 1; i < n; i++ {
		fmt.Fprintf(out, "%v\t%e\n", i, AlleleFrequencyProbability(i, n, alpha, binomCache, integralError))
	}

	err := out.Close()
	exception.PanicOnErr(err)
}

func PlotDerivedAscertainmentProbability(outFile string, n int, d int) {
	out := fileio.EasyCreate(outFile)
	fmt.Fprintf(out, "Frequency\tProbability\n")

	for i := 1; i < n; i++ {
		fmt.Fprintf(out, "%v\t%e\n", i, DerivedAscertainmentProbability(n, i, d))
	}

	err := out.Close()
	exception.PanicOnErr(err)
}

func PlotAncestralAscertainmentProbability(outFile string, n int, d int) {
	out := fileio.EasyCreate(outFile)
	fmt.Fprintf(out, "Frequency\tProbability\n")

	for i := 1; i < n; i++ {
		fmt.Fprintf(out, "%v\t%e\n", i, AncestralAscertainmentProbability(n, i, d))
	}

	err := out.Close()
	exception.PanicOnErr(err)
}

func PlotDerivedAscertainmentDenominator(outFile string, n int, d int, alpha float64, integralError float64) {
	allN := []int{n}
	binomCache := BuildBinomCache(allN)
	fCache := BuildFCache(n, alpha, binomCache, integralError)
	fCacheSum := GetFCacheSum(fCache)

	out := fileio.EasyCreate(outFile)
	fmt.Fprintf(out, "Frequency\tProbability\n")

	for i := 1; i < n; i++ {
		fmt.Fprintf(out, "%v\t%e\n", i, DerivedAscertainmentDenominator(fCache, fCacheSum, d))
	}

	err := out.Close()
	exception.PanicOnErr(err)
}

func PlotAncestralAscertainmentDenominator(outFile string, n int, d int, alpha float64, integralError float64) {
	allN := []int{n}
	binomCache := BuildBinomCache(allN)
	fCache := BuildFCache(n, alpha, binomCache, integralError)
	fCacheSum := GetFCacheSum(fCache)

	out := fileio.EasyCreate(outFile)
	fmt.Fprintf(out, "Frequency\tProbability\n")

	for i := 1; i < n; i++ {
		fmt.Fprintf(out, "%v\t%e\n", i, AncestralAscertainmentDenominator(fCache, fCacheSum, d))
	}

	err := out.Close()
	exception.PanicOnErr(err)
}
