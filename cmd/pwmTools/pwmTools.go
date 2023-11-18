package main

import (
	"flag"
	"fmt"
	"log"
)

func usage() {
	fmt.Printf(
		"pwmTools - a collection of tools for manipulating position matrices.\n" +
			"Usage:\n" +
			"pwmTools filter in.pwm out.pwm\n" +
			"OR\n" +
			"pwmTools info in.pwm out.txt\n" +
			"options:\n")
	flag.PrintDefaults()
}

func main() {
	flag.Usage = usage
	log.SetFlags(log.Ldate | log.Ltime | log.Lshortfile)
	flag.Parse()

	if len(flag.Args()) < 1 {
		flag.Usage()
		log.Fatalf("Error: user must specify a pwmTools subcommand.\n")
	}

	switch flag.Arg(0) {
	case "filter":
		parseFilterArgs()
	case "info":
		parseInfoArgs()
	default:
		flag.Usage()
		log.Fatalf("Error: unrecognized subcommand: %v.\n", flag.Arg(0))
	}
}
