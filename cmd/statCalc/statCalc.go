package main

import (
	"flag"
	"fmt"
	"github.com/vertgenlab/gonomics/common"
	"github.com/vertgenlab/gonomics/numbers"
	"log"
	"strings"
)

//kills the program if multiple options are selected.
func MultipleOptionErrorCheck(Normal *string, Binomial *string, Poisson *string, Beta *string, Gamma *string) {
	var count int = 0
	if *Normal != "" {
		count++
	}
	if *Binomial != "" {
		count++
	}
	if *Poisson != "" {
		count++
	}
	if *Beta != "" {
		count++
	}
	if *Gamma != "" {
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
			" -normal=mu,sigma. Defines a normal distribution with mean mu and standard deviation sigma. Ex Usage: -normal=0,1 1 or -normal=0,1 2 inf\n" +
			" -binomial=n,p. Defines a binomial distribution with n experiments and success probability p. Ex Usage: -binomial=10,0.5 3 or -binomial=10, 0.5 6 n\n" +
			" -poisson=lambda. Defines a poisson distribution with rate parameter lambda. Ex Usage: -poisson=4 4\n" +
			" -beta=alpha,beta. Defines a beta dsitribution with paramters alpha and beta. Ex Usage: -beta=5,5 0.2\n" +
			" -gamma=alpha,beta. Defines a gamma distribution with parameters alpha and beta. Ex Usage: -gamma=4,4 6\n" +
			"After defining a distribution, one float64 argument returns the function density at that value.\n" +
			"For discrete distributions, two arguments will evaluate the sum between two input values.\n" +
			"For the binomial distribution summation, the second argument can be set to n or N to evaluate the entire right tailed sum.\n" +
			"For continuous distributions, two arguments will evaluate an integral between the two input values with the defined distribution as the integrand.\n")
}

func main() {
	var Normal *string = flag.String("normal", "", "")
	var Binomial *string = flag.String("binomial", "", "")
	var Poisson *string = flag.String("poisson", "", "")
	var Beta *string = flag.String("beta", "", "")
	var Gamma *string = flag.String("gamma", "", "")

	flag.Usage = usage
	log.SetFlags(log.Ldate | log.Ltime | log.Lshortfile)
	flag.Parse()

	MultipleOptionErrorCheck(Normal, Binomial, Poisson, Beta, Gamma)

	if *Normal != "" {
		words := strings.Split(*Normal, ",")
		if len(words) != 2 {
			log.Fatalf("Error: A normal distribution is defined by two parameters. Received %v.\n", len(words))
		}
		mu := common.StringToFloat64(words[0])
		sigma := common.StringToFloat64(words[1])
		if len(flag.Args()) > 2 || len(flag.Args()) < 1 {
			flag.Usage()
			log.Fatalf("Error: expected one or two arguments, but got %d\n", len(flag.Args()))
		}
		if len(flag.Args()) == 1 {
			x := common.StringToFloat64(flag.Arg(0))
			fmt.Printf("%e\n", numbers.NormalDist(x, mu, sigma))
		} else if len(flag.Args()) == 2 {
			fmt.Printf("%e\n", numbers.NormalAdaptiveIntegral(flag.Arg(0), flag.Arg(1), mu, sigma))
		}
	} else if *Binomial != "" {
		words := strings.Split(*Binomial, ",")
		if len(words) != 2 {
			log.Fatalf("Error: a binomial distribution is defined by two parameters. Received %v.\n", len(words))
		}
		n := common.StringToInt(words[0])
		p := common.StringToFloat64(words[1])
		if len(flag.Args()) > 2 || len(flag.Args()) < 1 {
			flag.Usage()
			log.Fatalf("Error: expected one or two arguments, but got %d\n", len(flag.Args()))
		}
		if len(flag.Args()) == 1 {
			i := common.StringToInt(flag.Arg(0))
			fmt.Printf("%e\n", numbers.BinomialDist(n, i, p))
		} else if len(flag.Args()) == 2 {
			left := common.StringToInt(flag.Arg(0))
			if flag.Arg(1) == "N" || flag.Arg(1) == "n" {
				if left == 0 {
					fmt.Printf("%e\n", 1.00000)
				} else {
					fmt.Printf("%e\n", numbers.BinomialRightSummation(n, left, p))
				}
			} else if left == 0 {
				right := common.StringToInt(flag.Arg(1))
				fmt.Printf("%e\n", numbers.BinomialLeftSummation(n, right, p))
			} else {
				right := common.StringToInt(flag.Arg(1))
				fmt.Printf("%e\n", numbers.BinomialSum(left, right, n, p))
			}
		}
	} else if *Poisson != "" {
		lambda := common.StringToFloat64(*Poisson)
		if len(flag.Args()) > 2 || len(flag.Args()) < 1 {
			flag.Usage()
			log.Fatalf("Error: expected one or two arguments, but got %d\n", len(flag.Args()))
		}
		if len(flag.Args()) == 1 {
			k := common.StringToInt(flag.Arg(0))
			fmt.Printf("%e\n", numbers.PoissonDist(k, lambda))
		} else if len(flag.Args()) == 2 {
			if flag.Arg(1) == "INF" || flag.Arg(1) == "inf" || flag.Arg(1) == "Inf" {
				k := common.StringToInt(flag.Arg(0))
				fmt.Printf("%e\n", numbers.PoissonRightSummation(k, lambda))
			} else {
				left := common.StringToInt(flag.Arg(0))
				right := common.StringToInt(flag.Arg(1))
				fmt.Printf("%e\n", numbers.PoissonSum(left, right, lambda))
			}
		}
	} else if *Beta != "" {
		words := strings.Split(*Beta, ",")
		if len(words) != 2 {
			log.Fatalf("Error: a beta distribution is defined by two parameters. Received %v.\n", len(words))
		}
		alpha := common.StringToFloat64(words[0])
		beta := common.StringToFloat64(words[1])

		if len(flag.Args()) > 2 || len(flag.Args()) < 1 {
			flag.Usage()
			log.Fatalf("Error: expected one or two arguments, but got %d\n", len(flag.Args()))
		}
		if len(flag.Args()) == 1 {
			x := common.StringToFloat64(flag.Arg(0))
			fmt.Printf("%e\n", numbers.BetaDist(x, alpha, beta))
		} else if len(flag.Args()) == 2 {
			left := common.StringToFloat64(flag.Arg(0))
			right := common.StringToFloat64(flag.Arg(1))
			fmt.Printf("%e\n", numbers.BetaIntegral(left, right, alpha, beta))
		}
	} else if *Gamma != "" {
		words := strings.Split(*Gamma, ",")
		if len(words) != 2 {
			log.Fatalf("Error: a gamma distribution is defined by two parameters. Received %v.\n", len(words))
		}
		alpha := common.StringToFloat64(words[0])
		beta := common.StringToFloat64(words[1])
		if len(flag.Args()) > 2 || len(flag.Args()) < 1 {
			flag.Usage()
			log.Fatalf("Error: expected one or two arguments, but got %d\n", len(flag.Args()))
		}
		if len(flag.Args()) == 1 {
			x := common.StringToFloat64(flag.Arg(0))
			fmt.Printf("%e\n", numbers.GammaDist(x, alpha, beta))
		} else if len(flag.Args()) == 2 {
			left := common.StringToFloat64(flag.Arg(0))
			if flag.Arg(1) == "INF" || flag.Arg(1) == "inf" || flag.Arg(1) == "Inf" {
				fmt.Printf("%e\n", numbers.GammaRightIntegral(left, alpha, beta))
			} else {
				right := common.StringToFloat64(flag.Arg(1))
				fmt.Printf("%e\n", numbers.GammaIntegral(left, right, alpha, beta))
			}
		}
	} else {
		flag.Usage()
		log.Fatalf("Error: No distribution command specified.")
	}
}
