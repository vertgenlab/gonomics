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
			"\tpwmTools filter in.pwm out.pwm\n" +
			"\tOR\n" +
			"\tpwmTools format in.pwm out.pwm\n" +
			"\tOR\n" +
			"\tpwmTools info in.pwm out.txt\n" +
			"Enter a subcommand to view options.\n")
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
	case "format":
		parseFormatArgs()
	case "info":
		parseInfoArgs()
	default:
		flag.Usage()
		log.Fatalf("Error: unrecognized subcommand: %v.\n", flag.Arg(0))
	}
}
