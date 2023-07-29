package main

import "fmt"

func main() {
	fmt.Print(
		"This directory contains deprecated commands, preserved for backwards compatability.\n" +
			"Each deprecated command should point to an intended replacement program, and should\n" +
			"describe how the usage of a deprecated program can be transferred to usage in the new program.\n")
}
