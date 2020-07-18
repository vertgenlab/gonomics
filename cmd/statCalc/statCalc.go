package main

import (
	"log"
	"flag"
	"fmt"
	"github.com/vertgenlab/gonomics/numbers"
	"strconv"
)

func usage() {
	fmt.Print(
	"statCalc - Command line statistics calculator.\n" +
		"Usage:\t" +
		" statCalc [options] \n" +
		"options:\n" +
		" -normalDist x float64 mu float64 sigma float64\tCalculate the probability density for a given x value from a normal distribution with mean mu and standard deviation sigma.\n\tUsage Example: statCalc -normalDist 0 0 1\n" +
		" -binomialDist n int k int p float64\tCalculate the probability density for a given k value from a binomial distribution with n observations and success probability p. \n\tUsage Example: statCalc -binomialDist 10 2 0.5\n" +
		" -poissonDist k int lambda float64\tCalculate the probability density for a given k value from a poisson distribution with rate lambda. \n\tUsage Example: statCalc -poissonDist 5 5\n" + 
		" -betaDist x float64 alpha float64 beta float64\tCalculate the probability density for a given x value from a beta distribution with parameters alpha and beta. \n\tUsage Example: statCalc -betaDist 0.45 2.0 2.0\n" +
		" -gammaDist x float64 alpha float64 beta float64\tCalculate the probability density for a given x value from a gamma distribution with shape alpha and rate beta. \n\tUsage Example: statCalc -gammaDist 7.0 7.5 1.0\n" +
		" -poissonLeftSummation k int lambda float64\t Calculate the sum of probabilities in a Poisson distribution to the left of a particular k value, inclusive. \n\tUsage Example: statCalc -poissonLeftSummation 4 4\n" +
		" -poissonRightSummation k int lambda float64\t Calculate the sum of probabilities in a Poisson distribution to the right of a particular k value, inclusive. \n\tUsage Example: statCalc -poissonRightSummation 4 4\n" +
		" -binomialLeftSummation n int k int p float64\t Calculate the sum of probabilities in a binomial distribution to the left of a particular k value, inclusive. \n\tUsage Example: statCalc -binomialLeftSummation 10 2 0.5\n" +
		" -binomialRightSummation n int k int p float64\t Calculate the sum of probabilities in a binomial distribution to the right of a particular k value, inclusive. \n\tUsage Example: statCalc -binomialRightSummation 10 2 0.5\n" +
		" -normalLeftIntegral x float64 mu float64 sigma float64\t Calculate the integral of the normal distribution to the left of an input x value. \n\tUsage Example: statCalc -normalLeftIntegral -- -1.0 0.0 1.0\n" +
		" -normalRightIntegral x float64 mu float64 sigma float64\t Calculate the integral of the normal distribution to the right of an input x value. \n\tUsage Example: statCalc -normalRightIntegral -- -1.0 0.0 1.0\n" +
		" -betaLeftIntegral x float64 alpha float64 beta float64\t Calculate the integral fo the beta distribution to the left of an input x value. \n\tUsage Example: statCalc -betaLeftIntegral 0.2 1 1\n" +
		" -betaRightIntegral x float64 alpha float64 beta float64\t Calculate the integral fo the beta distribution to the right of an input x value. \n\tUsage Example: statCalc -betaRightIntegral 0.2 1 1\n" +
		" -gammaLeftIntegral x float64 alpha float64 beta float64\t Calculate the integral fo the gamma distribution to the left of an input x value. \n\tUsage Example: statCalc -gammaLeftIntegral 0.2 1 1\n" +
		" -gammaRightIntegral x float64 alpha float64 beta float64\t Calculate the integral fo the gamma distribution to the right of an input x value. \n\tUsage Example: statCalc -gammaRightIntegral 0.2 1 1\n")
}

func main() {
	var expectedNumArgs int
	var NormalDist *bool = flag.Bool("normalDist", false, "Calculate the probability density for a given x value from a normal distribution with mean mu and standard deviation sigma.")
	var BinomialDist *bool = flag.Bool("binomialDist", false, "Calculate the probability density for a given k value from a binomial distribution with n observations and success probability p.")
	var PoissonDist *bool = flag.Bool("poissonDist", false, "Calculate the probability density for a given k value from a poisson distribution with rate lambda.")
	var BetaDist *bool = flag.Bool("betaDist", false, "Calculate the probability density for a given x value from a beta distribution with parameters alpha and beta.")
	var GammaDist *bool = flag.Bool("gammaDist", false, "Calculate the probability density for a given x value from a gamma distribution with shape alpha and rate beta.")
	var PoissonLeftSummation *bool = flag.Bool("poissonLeftSummation", false, "Calculate the sum of probabilities in a Poisson distribution to the left of a particular k value, inclusive.")
	var PoissonRightSummation *bool = flag.Bool("poissonRightSummation", false,  "Calculate the sum of probabilities in a Poisson distribution to the right of a particular k value, inclusive.")
	var BinomialLeftSummation *bool = flag.Bool("binomialLeftSummation", false, "Calculate the sum of probabilities in a binomial distribution to the left of a particular k value, inclusive.")
	var BinomialRightSummation *bool = flag.Bool("binomialRightSummation", false, "Calculate the sum of probabilities in a binomial distribution to the right of a  particular k value, incluseive.")
	var NormalLeftIntegral *bool = flag.Bool("normalLeftIntegral", false, "Calculate the integral of the normal probability distribution to the left of an input x value.")
	var NormalRightIntegral *bool = flag.Bool("normalRightIntegral", false, "Calculate the integral of the normal probability distribution to the right of an input x value.")
	var BetaLeftIntegral *bool = flag.Bool("betaLeftIntegral", false, "Calculate the integral of the beta distribution to the left of an input x value.")
	var BetaRightIntegral *bool = flag.Bool("betaRightIntegral", false, "Calculate the integral of the beta distribution to the right of an input x value.")
	var GammaLeftIntegral *bool = flag.Bool("gammaLeftIntegral", false, "Calculate the integral of the gamma distribution to the left of an input x value.")
	var GammaRightIntegral *bool = flag.Bool("gammaRightIntegral", false, "Calculate the integral of the gamma distribution to the right of an input x value.")

	flag.Usage = usage
	log.SetFlags(log.Ldate | log.Ltime | log.Lshortfile)
	flag.Parse()
	
	if *NormalDist {
		expectedNumArgs = 3
		if len(flag.Args()) != expectedNumArgs {
			flag.Usage()
			log.Fatalf("Error: expecting %d arguments, but got %d\n",
				expectedNumArgs, len(flag.Args()))
		}
		x, _ := strconv.ParseFloat(flag.Arg(0), 64)
		mu, _ := strconv.ParseFloat(flag.Arg(1), 64)
		sigma, _ := strconv.ParseFloat(flag.Arg(2), 64)
		fmt.Printf("%e\n", numbers.NormalDist(x, mu, sigma))
	} else if *BinomialDist {
		expectedNumArgs = 3
		if len(flag.Args()) != expectedNumArgs {
			flag.Usage()
			log.Fatalf("Error: expecting %d arguments, but got %d\n",
				expectedNumArgs, len(flag.Args()))
		}
		n, _ := strconv.Atoi(flag.Arg(0))
		i, _ := strconv.Atoi(flag.Arg(1))
		p, _ := strconv.ParseFloat(flag.Arg(2), 64)
		fmt.Printf("%e\n", numbers.BinomialDist(n, i, p))
	} else if *PoissonDist {
		expectedNumArgs = 2
		if len(flag.Args()) != expectedNumArgs {
			flag.Usage()
			log.Fatalf("Error: expecting %d arguments, but got %d\n",
				expectedNumArgs, len(flag.Args()))
		}
		k, _ := strconv.Atoi(flag.Arg(0))
		lambda, _ := strconv.ParseFloat(flag.Arg(1), 64)
		fmt.Printf("%e\n", numbers.PoissonDist(k, lambda))
	} else if *BetaDist {
		expectedNumArgs = 3
		if len(flag.Args()) != expectedNumArgs {
			flag.Usage()
			log.Fatalf("Error: expecting %d arguments, but got %d\n",
				expectedNumArgs, len(flag.Args()))
		}
		x, _ := strconv.ParseFloat(flag.Arg(0), 64)
		alpha, _ := strconv.ParseFloat(flag.Arg(1), 64)
		beta, _ := strconv.ParseFloat(flag.Arg(2), 64)
		fmt.Printf("%e\n", numbers.BetaDist(x, alpha, beta))
	} else if *GammaDist {
		expectedNumArgs = 3
		if len(flag.Args()) != expectedNumArgs {
			flag.Usage()
			log.Fatalf("Error: expecting %d arguments, but got %d\n",
				expectedNumArgs, len(flag.Args()))
		}
		x, _ := strconv.ParseFloat(flag.Arg(0), 64)
		alpha, _ := strconv.ParseFloat(flag.Arg(1), 64)
		beta, _ := strconv.ParseFloat(flag.Arg(2), 64)
		fmt.Printf("%e\n", numbers.GammaDist(x, alpha, beta))
	} else if *PoissonLeftSummation {
		expectedNumArgs = 2
		if len(flag.Args()) != expectedNumArgs {
			flag.Usage()
			log.Fatalf("Error: expecting %d arguments, but got %d\n",
				expectedNumArgs, len(flag.Args()))
		}
		k, _ := strconv.Atoi(flag.Arg(0))
		lambda, _ := strconv.ParseFloat(flag.Arg(1), 64)
		fmt.Printf("%e\n", numbers.PoissonLeftSummation(k, lambda))
	} else if *PoissonRightSummation {
		expectedNumArgs = 2
		if len(flag.Args()) != expectedNumArgs {
			flag.Usage()
			log.Fatalf("Error: expecting %d arguments, but got %d\n",
				expectedNumArgs, len(flag.Args()))
		}
		k, _ := strconv.Atoi(flag.Arg(0))
		lambda, _ := strconv.ParseFloat(flag.Arg(1), 64)
		fmt.Printf("%e\n", numbers.PoissonRightSummation(k, lambda))
	} else if *BinomialLeftSummation {
		expectedNumArgs = 3
		if len(flag.Args()) != expectedNumArgs {
			flag.Usage()
			log.Fatalf("Error: expecting %d arguments, but got %d\n",
				expectedNumArgs, len(flag.Args()))
		}
		n, _ := strconv.Atoi(flag.Arg(0))
		i, _ := strconv.Atoi(flag.Arg(1))
		p, _ := strconv.ParseFloat(flag.Arg(2), 64)
		fmt.Printf("%e\n", numbers.BinomialLeftSummation(n, i, p))
	} else if *BinomialRightSummation {
		expectedNumArgs = 3
		if len(flag.Args()) != expectedNumArgs {
			flag.Usage()
			log.Fatalf("Error: expecting %d arguments, but got %d\n",
				expectedNumArgs, len(flag.Args()))
		}
		n, _ := strconv.Atoi(flag.Arg(0))
		i, _ := strconv.Atoi(flag.Arg(1))
		p, _ := strconv.ParseFloat(flag.Arg(2), 64)
		fmt.Printf("%e\n", numbers.BinomialRightSummation(n, i, p))
	} else if *NormalLeftIntegral {
		//integrates to 200 standard deviations above the mean
		expectedNumArgs = 3
		x, _ := strconv.ParseFloat(flag.Arg(0), 64)
		mu, _ := strconv.ParseFloat(flag.Arg(1), 64)
		sigma, _ := strconv.ParseFloat(flag.Arg(2), 64)
		f := numbers.NormalClosure(mu, sigma)
		leftBound := mu - 200.0 * sigma
		fmt.Printf("%e\n", numbers.DefiniteIntegral(f, leftBound, x))
	} else if *NormalRightIntegral {
		//integrates to 200 standard deviations above the mean
		expectedNumArgs = 3
		x, _ := strconv.ParseFloat(flag.Arg(0), 64)
		mu, _ := strconv.ParseFloat(flag.Arg(1), 64)
		sigma, _ := strconv.ParseFloat(flag.Arg(2), 64)
		f := numbers.NormalClosure(mu, sigma)
		rightBound := mu + 200.0 * sigma
		fmt.Printf("%e\n", numbers.DefiniteIntegral(f, x, rightBound))
	} else if *BetaLeftIntegral {
		expectedNumArgs = 3
		if len(flag.Args()) != expectedNumArgs {
			flag.Usage()
			log.Fatalf("Error: expecting %d arguments, but got %d\n",
				expectedNumArgs, len(flag.Args()))
		}
		x, _ := strconv.ParseFloat(flag.Arg(0), 64)
		alpha, _ := strconv.ParseFloat(flag.Arg(1), 64)
		beta, _ := strconv.ParseFloat(flag.Arg(2), 64)
		f := numbers.BetaClosure(alpha, beta)
		fmt.Printf("%e\n", numbers.DefiniteIntegral(f, 0, x))
	} else if *BetaRightIntegral {
		expectedNumArgs = 3
		if len(flag.Args()) != expectedNumArgs {
			flag.Usage()
			log.Fatalf("Error: expecting %d arguments, but got %d\n",
				expectedNumArgs, len(flag.Args()))
		}
		x, _ := strconv.ParseFloat(flag.Arg(0), 64)
		alpha, _ := strconv.ParseFloat(flag.Arg(1), 64)
		beta, _ := strconv.ParseFloat(flag.Arg(2), 64)
		f := numbers.BetaClosure(alpha, beta)
		fmt.Printf("%e\n", numbers.DefiniteIntegral(f, x, 1))
	} else if *GammaLeftIntegral {
		expectedNumArgs = 3
		if len(flag.Args()) != expectedNumArgs {
			flag.Usage()
			log.Fatalf("Error: expecting %d arguments, but got %d\n",
				expectedNumArgs, len(flag.Args()))
		}
		x, _ := strconv.ParseFloat(flag.Arg(0), 64)
		alpha, _ := strconv.ParseFloat(flag.Arg(1), 64)
		beta, _ := strconv.ParseFloat(flag.Arg(2), 64)
		f := numbers.GammaClosure(alpha, beta)
		fmt.Printf("%e\n", numbers.DefiniteIntegral(f, 0, x))
	} else if *GammaRightIntegral { 
		expectedNumArgs = 3
		if len(flag.Args()) != expectedNumArgs {
			flag.Usage()
			log.Fatalf("Error: expecting %d arguments, but got %d\n",
				expectedNumArgs, len(flag.Args()))
		}
		x, _ := strconv.ParseFloat(flag.Arg(0), 64)
		alpha, _ := strconv.ParseFloat(flag.Arg(1), 64)
		beta, _ := strconv.ParseFloat(flag.Arg(2), 64)
		f := numbers.GammaClosure(alpha, beta)
		fmt.Printf("%e\n", 1 - numbers.DefiniteIntegral(f, 0, x))
	} else {
		flag.Usage()
		log.Fatalf("Error: No distribution command specified.")
	}
}