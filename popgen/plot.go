package popgen

import (
	"fmt"
	"log"

	"github.com/vertgenlab/gonomics/exception"
	"github.com/vertgenlab/gonomics/fileio"
)

// PlotAfsF writes the Allele Frequency F function (AfsSampleDensity) to an output file for downstream visualization.
func PlotAfsF(alpha float64, n int, outFile string, integralError float64) {
	var err error
	out := fileio.EasyCreate(outFile)
	allN := []int{n}
	binomCache := BuildBinomCache(allN)

	_, err = fmt.Fprintf(out, "Frequency\tF\n")
	exception.PanicOnErr(err)

	for i := 1; i < n; i++ {
		_, err = fmt.Fprintf(out, "%v\t%e\n", i, AfsSampleDensity(n, i, alpha, binomCache, integralError))
		exception.PanicOnErr(err)
	}
	err = out.Close()
	exception.PanicOnErr(err)
}

// PlotAfsPmf writes the allele frequency probability mass function (AlleleFrequencyProbability) to an output file for downstream visualization.
// Derived or ancestral flag enable visualization of the AfsPmf with the divergence-based ascertainment correction.
func PlotAfsPmf(alpha float64, n int, outFile string, integralError float64, derived bool, ancestral bool) {
	if derived && ancestral {
		log.Fatalf("PlotAfsPmf cannot be used with both the derived and ancestral correction simultaneously.")
	}
	var err error
	out := fileio.EasyCreate(outFile)

	allN := []int{n}
	binomCache := BuildBinomCache(allN)

	//a short header is written to the file
	_, err = fmt.Fprintf(out, "Frequency\tProbability\n")
	exception.PanicOnErr(err)
	//for each possible allele frequency for a segregating site, we calculate the probability to the file
	for i := 1; i < n; i++ {
		if ancestral {
			_, err = fmt.Fprintf(out, "%v\t%e\n", i, AlleleFrequencyProbabilityAncestralAscertainment(alpha, i, n, 1, binomCache, integralError))
		} else if derived {
			_, err = fmt.Fprintf(out, "%v\t%e\n", i, AlleleFrequencyProbabilityDerivedAscertainment(alpha, i, n, 1, binomCache, integralError))
		} else {
			_, err = fmt.Fprintf(out, "%v\t%e\n", i, AlleleFrequencyProbability(i, n, alpha, binomCache, integralError))
		}
		exception.PanicOnErr(err)
	}

	err = out.Close()
	exception.PanicOnErr(err)
}

// PlotAfsLikelihood plots the likelihood function from input Afs data. Specifies a left and right bound for plotting and a number of points to be plotted.
// Optiosn available for applying the divergence based ascertainment bias correction, changing the ascertainment subset size, and changing the integral error.
// This program will fatal if it is asked to plot the likelihood when alpha is equal to exactly 0.
func PlotAfsLikelihood(afs Afs, outFile string, leftBound float64, rightBound float64, numPoints int, integralError float64, divergenceAscertainment bool, D int) {
	var err error
	var alpha float64
	out := fileio.EasyCreate(outFile)

	allN := findAllN(afs)
	binomMap := BuildBinomCache(allN)

	_, err = fmt.Fprintf(out, "Alpha\tLikelihood\n")
	exception.PanicOnErr(err)

	for i := 0; i <= numPoints; i++ {
		alpha = leftBound + (float64(i)/float64(numPoints))*(rightBound-leftBound)
		if divergenceAscertainment {
			_, err = fmt.Fprintf(out, "%e\t%e\n", alpha, AfsDivergenceAscertainmentFixedAlpha(afs, alpha, binomMap, D, integralError))
			exception.PanicOnErr(err)
		} else {
			_, err = fmt.Fprintf(out, "%e\t%e\n", alpha, AfsLikelihoodFixedAlpha(afs, alpha, binomMap, integralError))
			exception.PanicOnErr(err)
		}
	}

	err = out.Close()
	exception.PanicOnErr(err)
}

func PlotDerivedAscertainmentProbability(outFile string, n int, d int) {
	var err error
	out := fileio.EasyCreate(outFile)
	_, err = fmt.Fprintf(out, "Frequency\tProbability\n")
	exception.PanicOnErr(err)

	for i := 1; i < n; i++ {
		_, err = fmt.Fprintf(out, "%v\t%e\n", i, DerivedAscertainmentProbability(n, i, d))
		exception.PanicOnErr(err)
	}

	err = out.Close()
	exception.PanicOnErr(err)
}

func PlotAncestralAscertainmentProbability(outFile string, n int, d int) {
	var err error
	out := fileio.EasyCreate(outFile)
	_, err = fmt.Fprintf(out, "Frequency\tProbability\n")
	exception.PanicOnErr(err)

	for i := 1; i < n; i++ {
		_, err = fmt.Fprintf(out, "%v\t%e\n", i, AncestralAscertainmentProbability(n, i, d))
		exception.PanicOnErr(err)
	}

	err = out.Close()
	exception.PanicOnErr(err)
}

func PlotDerivedAscertainmentDenominator(outFile string, n int, d int, alpha float64, integralError float64) {
	var err error
	allN := []int{n}
	binomCache := BuildBinomCache(allN)
	fCache := BuildFCache(n, alpha, binomCache, integralError)
	fCacheSum := GetFCacheSum(fCache)

	out := fileio.EasyCreate(outFile)
	_, err = fmt.Fprintf(out, "Frequency\tProbability\n")
	exception.PanicOnErr(err)

	for i := 1; i < n; i++ {
		_, err = fmt.Fprintf(out, "%v\t%e\n", i, DerivedAscertainmentDenominator(fCache, fCacheSum, d))
		exception.PanicOnErr(err)
	}

	err = out.Close()
	exception.PanicOnErr(err)
}

func PlotAncestralAscertainmentDenominator(outFile string, n int, d int, alpha float64, integralError float64) {
	var err error
	allN := []int{n}
	binomCache := BuildBinomCache(allN)
	fCache := BuildFCache(n, alpha, binomCache, integralError)
	fCacheSum := GetFCacheSum(fCache)

	out := fileio.EasyCreate(outFile)
	_, err = fmt.Fprintf(out, "Frequency\tProbability\n")
	exception.PanicOnErr(err)

	for i := 1; i < n; i++ {
		_, err = fmt.Fprintf(out, "%v\t%e\n", i, AncestralAscertainmentDenominator(fCache, fCacheSum, d))
		exception.PanicOnErr(err)
	}

	err = out.Close()
	exception.PanicOnErr(err)
}
