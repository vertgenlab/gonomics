package popgen

import (
	"fmt"
	"github.com/vertgenlab/gonomics/exception"
	"github.com/vertgenlab/gonomics/fileio"
)

func PlotAfsF(alpha float64, n int, outFile string) {
	out := fileio.EasyCreate(outFile)
	allN := []int{n}
	binomCache := BuildBinomCache(allN)

	fmt.Fprintf(out, "Frequency\tF\n")

	for i := 1; i < n; i++ {
		fmt.Fprintf(out, "%v\t%e\n", i, AFSSampleDensity(n, i, alpha, binomCache))
	}
	err := out.Close()
	exception.PanicOnErr(err)
}

func PlotAfsPmf(alpha float64, n int, outFile string) {
	out := fileio.EasyCreate(outFile)

	allN := []int{n}
	binomCache := BuildBinomCache(allN)
	//DEBUG:fmt.Printf("Len BinomCache:%v. Entry at index n is of length: %v.\n", len(binomCache), len(binomCache[n]))

	//a short header is written to the file
	fmt.Fprintf(out, "Frequency\tProbability\n")
	//for each possible allele frequency for a segregating site, we calculate the probability to the file
	for i := 1; i < n; i++ {
		fmt.Fprintf(out, "%v\t%e\n", i, AlleleFrequencyProbability(i, n, alpha, binomCache))
	}

	err := out.Close()
	exception.PanicOnErr(err)
}
