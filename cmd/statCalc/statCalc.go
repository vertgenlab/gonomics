// Command Group: "Statistics & Population Genetics"

package main

import (
	"flag"
	"fmt"
	"log"
	"math/rand"
	"strings"

	"github.com/vertgenlab/gonomics/common"
	"github.com/vertgenlab/gonomics/exception"
	"github.com/vertgenlab/gonomics/fileio"
	"github.com/vertgenlab/gonomics/numbers"
	"github.com/vertgenlab/gonomics/popgen"
)

func statCalc(s Settings) {
	MultipleOptionErrorCheck(s.Normal, s.Binomial, s.Poisson, s.Beta, s.Gamma, s.SampleAfs, s.SampleBeta, s.SampleGamma, s.SampleNormal)
	rand.Seed(s.SetSeed)
	var err error
	out := fileio.EasyCreate(s.OutFile)
	if s.Normal != "" {
		words := strings.Split(s.Normal, ",")
		if len(words) != 2 {
			log.Fatalf("Error: A normal distribution is defined by two parameters. Received %v.\n", len(words))
		}
		mu := common.StringToFloat64(words[0])
		sigma := common.StringToFloat64(words[1])
		if len(s.Args) > 2 || len(s.Args) < 1 {
			flag.Usage()
			log.Fatalf("Error: expected one or two arguments, but got %d\n", len(s.Args))
		}
		if len(s.Args) == 1 {
			x := common.StringToFloat64(s.Args[0])
			_, err = fmt.Fprintf(out, "%e\n", numbers.NormalDist(x, mu, sigma))
			exception.PanicOnErr(err)
		} else if len(s.Args) == 2 {
			_, err = fmt.Fprintf(out, "%e\n", numbers.NormalAdaptiveIntegral(s.Args[0], s.Args[1], mu, sigma))
			exception.PanicOnErr(err)
		}
	} else if s.Binomial != "" {
		words := strings.Split(s.Binomial, ",")
		if len(words) != 2 {
			log.Fatalf("Error: a binomial distribution is defined by two parameters. Received %v.\n", len(words))
		}
		n := common.StringToInt(words[0])
		p := common.StringToFloat64(words[1])
		if len(s.Args) > 2 || len(s.Args) < 1 {
			flag.Usage()
			log.Fatalf("Error: expected one or two arguments, but got %d\n", len(s.Args))
		}
		if len(s.Args) == 1 {
			i := common.StringToInt(s.Args[0])
			answer, _ := numbers.BinomialDist(n, i, p)
			_, err = fmt.Fprintf(out, "%e\n", answer)
			exception.PanicOnErr(err)
		} else if len(s.Args) == 2 {
			left := common.StringToInt(s.Args[0])
			if s.Args[1] == "N" || s.Args[1] == "n" {
				if left == 0 {
					_, err = fmt.Fprintf(out, "%e\n", 1.00000)
					exception.PanicOnErr(err)
				} else {
					_, err = fmt.Fprintf(out, "%e\n", numbers.BinomialRightSummation(n, left, p))
					exception.PanicOnErr(err)
				}
			} else if left == 0 {
				right := common.StringToInt(s.Args[1])
				_, err = fmt.Fprintf(out, "%e\n", numbers.BinomialLeftSummation(n, right, p))
				exception.PanicOnErr(err)
			} else {
				right := common.StringToInt(s.Args[1])
				_, err = fmt.Fprintf(out, "%e\n", numbers.BinomialSum(left, right, n, p))
				exception.PanicOnErr(err)
			}
		}
	} else if s.Poisson != "" {
		lambda := common.StringToFloat64(s.Poisson)
		if len(s.Args) > 2 || len(s.Args) < 1 {
			flag.Usage()
			log.Fatalf("Error: expected one or two arguments, but got %d\n", len(flag.Args()))
		}
		if len(s.Args) == 1 {
			k := common.StringToInt(s.Args[0])
			_, err = fmt.Fprintf(out, "%e\n", numbers.PoissonDist(k, lambda))
			exception.PanicOnErr(err)
		} else if len(s.Args) == 2 {
			if s.Args[1] == "INF" || s.Args[1] == "inf" || s.Args[1] == "Inf" {
				k := common.StringToInt(s.Args[0])
				_, err = fmt.Fprintf(out, "%e\n", numbers.PoissonRightSummation(k, lambda))
				exception.PanicOnErr(err)
			} else {
				left := common.StringToInt(s.Args[0])
				right := common.StringToInt(s.Args[1])
				_, err = fmt.Fprintf(out, "%e\n", numbers.PoissonSum(left, right, lambda))
				exception.PanicOnErr(err)
			}
		}
	} else if s.Beta != "" {
		words := strings.Split(s.Beta, ",")
		if len(words) != 2 {
			log.Fatalf("Error: a beta distribution is defined by two parameters. Received %v.\n", len(words))
		}
		alpha := common.StringToFloat64(words[0])
		beta := common.StringToFloat64(words[1])

		if len(s.Args) > 2 || len(s.Args) < 1 {
			flag.Usage()
			log.Fatalf("Error: expected one or two arguments, but got %d\n", len(flag.Args()))
		}
		if len(s.Args) == 1 {
			x := common.StringToFloat64(s.Args[0])
			_, err = fmt.Fprintf(out, "%e\n", numbers.BetaDist(x, alpha, beta))
			exception.PanicOnErr(err)
		} else if len(s.Args) == 2 {
			left := common.StringToFloat64(s.Args[0])
			right := common.StringToFloat64(s.Args[1])
			_, err = fmt.Fprintf(out, "%e\n", numbers.BetaIntegral(left, right, alpha, beta))
			exception.PanicOnErr(err)
		}
	} else if s.Gamma != "" {
		words := strings.Split(s.Gamma, ",")
		if len(words) != 2 {
			log.Fatalf("Error: a gamma distribution is defined by two parameters. Received %v.\n", len(words))
		}
		alpha := common.StringToFloat64(words[0])
		beta := common.StringToFloat64(words[1])
		if len(s.Args) > 2 || len(s.Args) < 1 {
			flag.Usage()
			log.Fatalf("Error: expected one or two arguments, but got %d\n", len(s.Args))
		}
		if len(s.Args) == 1 {
			x := common.StringToFloat64(s.Args[0])
			_, err = fmt.Fprintf(out, "%e\n", numbers.GammaDist(x, alpha, beta))
			exception.PanicOnErr(err)
		} else if len(s.Args) == 2 {
			left := common.StringToFloat64(s.Args[0])
			if s.Args[1] == "INF" || s.Args[1] == "inf" || s.Args[1] == "Inf" {
				_, err = fmt.Fprintf(out, "%e\n", numbers.GammaRightIntegral(left, alpha, beta))
				exception.PanicOnErr(err)
			} else {
				right := common.StringToFloat64(s.Args[1])
				_, err = fmt.Fprintf(out, "%e\n", numbers.GammaIntegral(left, right, alpha, beta))
				exception.PanicOnErr(err)
			}
		}
	} else if s.SampleAfs != "" {
		words := strings.Split(s.SampleAfs, ",")
		if len(words) != 6 {
			log.Fatalf("Error: sampleAFS expected six parameters, received: %v.\n", len(words))
		}
		alpha := common.StringToFloat64(words[0])
		numSamples := common.StringToInt(words[1])
		maxSampleDepth := common.StringToInt(words[2])
		bins := common.StringToInt(words[3])
		xLeft := common.StringToFloat64(words[4])
		xRight := common.StringToFloat64(words[5])
		answer := popgen.StationaritySampler(alpha, numSamples, maxSampleDepth, bins, xLeft, xRight)
		for i := 0; i < len(answer); i++ {
			_, err = fmt.Fprintf(out, "%e\n", answer[i])
			exception.PanicOnErr(err)
		}
	} else if s.SampleBeta != "" {
		words := strings.Split(s.SampleBeta, ",")
		if len(words) != 3 {
			log.Fatalf("Error: sampleBeta expected four parameters, received: %v.\n", len(words))
		}
		alpha := common.StringToFloat64(words[0])
		beta := common.StringToFloat64(words[1])
		numSamples := common.StringToInt(words[2])
		sampler := numbers.BetaSampler(alpha, beta)
		var current float64
		for i := 0; i < numSamples; i++ {
			current, _ = sampler()
			_, err = fmt.Fprintf(out, "%e\n", current)
			exception.PanicOnErr(err)
		}
	} else if s.SampleGamma != "" {
		words := strings.Split(s.SampleGamma, ",")
		if len(words) != 3 {
			log.Fatalf("Error: sampleGamma expected four parameters, received: %v.\n", len(words))
		}
		alpha := common.StringToFloat64(words[0])
		beta := common.StringToFloat64(words[1])
		numSamples := common.StringToInt(words[2])
		sampler := numbers.GammaSampler(alpha, beta)
		var current float64
		for i := 0; i < numSamples; i++ {
			current, _ = sampler()
			_, err = fmt.Fprintf(out, "%e\n", current)
			exception.PanicOnErr(err)
		}
	} else if s.SampleNormal != "" {
		words := strings.Split(s.SampleNormal, ",")
		if len(words) != 3 {
			log.Fatalf("Error: sampleNormal expected four parameters, received: %v.\n", len(words))
		}
		mu := common.StringToFloat64(words[0])
		sigma := common.StringToFloat64(words[1])
		numSamples := common.StringToInt(words[2])
		for i := 0; i < numSamples; i++ {
			_, err = fmt.Fprintf(out, "%e\n", numbers.SampleInverseNormal(mu, sigma))
			exception.PanicOnErr(err)
		}
	} else {
		flag.Usage()
		log.Fatalf("Error: No distribution command specified.")
	}
	err = out.Close()
	exception.PanicOnErr(err)
}

// kills the program if multiple options are selected.
func MultipleOptionErrorCheck(Normal string, Binomial string, Poisson string, Beta string, Gamma string, SampleAfs string, SampleBeta string, SampleGamma string, SampleNormal string) {
	var count int = 0
	if Normal != "" {
		count++
	}
	if Binomial != "" {
		count++
	}
	if Poisson != "" {
		count++
	}
	if Beta != "" {
		count++
	}
	if Gamma != "" {
		count++
	}
	if SampleAfs != "" {
		count++
	}
	if SampleBeta != "" {
		count++
	}
	if SampleGamma != "" {
		count++
	}
	if SampleNormal != "" {
		count++
	}
	if count > 1 {
		log.Fatalf("Error: Multiple distribution arguments selected.")
	}
}

func usage() {
	fmt.Print(
		"statCalc - Command line statistics calculator.\n" +
			"Usage:\t" +
			" statCalc [options] \n" +
			"options:\n" +
			"-outFile. Redirect the output to a file. Default to stdout.\n\n" +
			"Enter one of the following options:\n" +
			" -normal=mu,sigma. Defines a normal distribution with mean mu and standard deviation sigma. Ex Usage: -normal=0,1 1 or -normal=0,1 2 inf\n" +
			" -binomial=n,p. Defines a binomial distribution with n experiments and success probability p. Ex Usage: -binomial=10,0.5 3 or -binomial=10, 0.5 6 n\n" +
			" -poisson=lambda. Defines a poisson distribution with rate parameter lambda. Ex Usage: -poisson=4 4\n" +
			" -beta=alpha,beta. Defines a beta dsitribution with parameters alpha and beta. Ex Usage: -beta=5,5 0.2\n" +
			" -gamma=alpha,beta. Defines a gamma distribution with parameters alpha and beta. Ex Usage: -gamma=4,4 6\n" +
			" -sampleAfs=alpha,numSamples,maxSampleDepth,bins,xLeft,xRight. Provides a list of values sampled from an allele frequency spectrum with selection parameter alpha.\n" +
			"\tsampleAFS will return numSamples many values between xLeft and xRight. Bins and maxSampleDepth are performance and accuracy options, suggested values are 1000 and 1000, respectively.\n" +
			"\tAfter defining a distribution, one float64 argument returns the function density at that value. Ex usage: -sampleAfs=0.02,200,1000,1000,0.001,0.999\n" +
			" -sampleBeta=alpha,beta,numSamples. Provides a list of values sampled from the beta distribution with a selected alpha and beta parameter.\n" +
			" -sampleGamma=alpha,beta,numSample. Provides a list of values sampled from the gamma distribution with a selected alpha and beta parameter.\n" +
			" -sampleNormal=mu,sigma,numSamples. Provides a list of values sampled from the normal distribution with a selected mu and sigma parameter.\n\n" +
			"For discrete distributions, two arguments will evaluate the sum between two input values.\n" +
			"For the binomial distribution summation, the second argument can be set to n or N to evaluate the entire right tailed sum.\n" +
			"For continuous distributions, two arguments will evaluate an integral between the two input values with the defined distribution as the integrand.\n")
}

type Settings struct {
	Args         []string
	Normal       string
	Binomial     string
	Poisson      string
	Beta         string
	Gamma        string
	SampleAfs    string
	SampleBeta   string
	SampleGamma  string
	SampleNormal string
	SetSeed      int64
	OutFile      string
}

func main() {
	var Normal *string = flag.String("normal", "", "")
	var Binomial *string = flag.String("binomial", "", "")
	var Poisson *string = flag.String("poisson", "", "")
	var Beta *string = flag.String("beta", "", "")
	var Gamma *string = flag.String("gamma", "", "")
	var SampleAfs *string = flag.String("sampleAfs", "", "")
	var SampleBeta *string = flag.String("sampleBeta", "", "")
	var SampleGamma *string = flag.String("sampleGamma", "", "")
	var SampleNormal *string = flag.String("sampleNormal", "", "")
	var setSeed *int64 = flag.Int64("setSeed", -1, "Use a specific seed for the RNG.")
	var outFile *string = flag.String("outFile", "stdout", "Redirect the output to a file. Default to stdout.")

	flag.Usage = usage
	//log.SetFlags(log.Ldate | log.Ltime | log.Lshortfile)
	flag.Parse()

	s := Settings{
		Args:         flag.Args(),
		Normal:       *Normal,
		Binomial:     *Binomial,
		Poisson:      *Poisson,
		Beta:         *Beta,
		Gamma:        *Gamma,
		SampleAfs:    *SampleAfs,
		SampleBeta:   *SampleBeta,
		SampleGamma:  *SampleGamma,
		SampleNormal: *SampleNormal,
		SetSeed:      *setSeed,
		OutFile:      *outFile,
	}

	statCalc(s)
}
