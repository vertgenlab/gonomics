package main

import (
	"flag"
	"fmt"
	"github.com/vertgenlab/gonomics/exception"
	"io/ioutil"
	"log"
	"os"
	"os/exec"
	"sort"
	"strings"
)

func usage() {
	fmt.Print(
		"gonomics - A collection of builtin tools that use the gonomics core library.\n\n" +
			"Usage: gonomics <command> [options]\n\n")
	printCmdList()
}

func getBin() (path string, binExists map[string]bool) {
	path = os.Getenv("GOBIN")

	if path == "" { // if GOBIN is not set
		//fmt.Println("WARNING: GOBIN is not set, defaulting to $GOPATH/bin")
		path = os.Getenv("GOPATH") + "/bin" // default to $GOPATH/bin
	}

	dir, err := ioutil.ReadDir(path)
	if err != nil {
		log.Printf("ERROR: could not access %s\n", path)
		log.Fatal(err)
	}

	binExists = make(map[string]bool)

	for _, file := range dir {
		binExists[file.Name()] = true
	}

	if !binExists["gonomics"] { // check if gonomics is installed by checking for this binary
		log.Fatalf("ERROR: gonomics binary were not found in %s\n "+
			"Run 'go install ./...' in gonomics directory to compile binary", path)
	}

	return path, binExists
}

func getGonomicsCmds() map[string]bool {
	expectedPath := os.Getenv("GOPATH") + "/src/github.com/vertgenlab/gonomics/cmd"
	cmds, err := ioutil.ReadDir(expectedPath)
	if err != nil {
		log.Printf("ERROR: could not find gonomics cmd folder in expected path: %s\n", expectedPath)
		log.Fatal(err)
	}

	funcNames := make(map[string]bool)

	for _, x := range cmds {
		funcNames[x.Name()] = true
	}

	if !funcNames["gonomics"] { // check that the directory for this cmd was found
		log.Fatalf("ERROR: gonomics cmd folder was not found in %s\n", expectedPath)
	}

	return funcNames
}

func printCmdList() {
	binPath, _ := getBin()
	cmdMap := getGonomicsCmds()

	cmds := make([]string, 0, len(cmdMap))
	for key := range cmdMap {
		cmds = append(cmds, key)
	}
	sort.Slice(cmds, func(i, j int) bool { return cmds[i] < cmds[j] })

	var endFirstLineIdx int
	var rawOutput []byte

	fmt.Println("Commands:")
	for _, cmdName := range cmds {
		cmd := exec.Command(binPath + "/" + cmdName)
		rawOutput, _ = cmd.Output()
		endFirstLineIdx = strings.Index(string(rawOutput), "\n")

		switch {
		case endFirstLineIdx > 0: // cmd starts with summary line
			fmt.Printf("     %s\n", string(rawOutput[:endFirstLineIdx]))

		case len(rawOutput) == 0: // cmd has no usage statement
			fmt.Printf("     %s\n", cmdName)

		case endFirstLineIdx == -1 && len(rawOutput) > 0: // has output, but no newline
			fmt.Printf("     %s\n", rawOutput)

		default: // if all else fails, print cmd name
			fmt.Printf("     %s\n", cmdName)
		}
	}
}

func main() {
	flag.Usage = usage
	log.SetFlags(log.Ldate | log.Ltime | log.Lshortfile)
	flag.Parse()

	if len(flag.Args()) == 0 {
		flag.Usage()
		return
	}

	cmdCalled := flag.Arg(0) // command called by the user (e.g. 'gonomics faFormat')

	binPath, binExists := getBin()
	isGonomicsFunc := getGonomicsCmds()

	switch { // Error checks. verbose for clarity
	case isGonomicsFunc[cmdCalled] && binExists[cmdCalled]:
		// proceed in main function

	case isGonomicsFunc[cmdCalled] && !binExists[cmdCalled]:
		log.Fatalf("ERROR: binary for %s was not found. "+
			"Please run 'go install ./...' in the gonomics directory.", cmdCalled)

	case !isGonomicsFunc[cmdCalled] && binExists[cmdCalled]:
		log.Fatalf("ERROR: %s is not a gonomics function", cmdCalled)

	default:
		log.Fatalf("ERROR: %s does not exist, and is not a gonomics function", cmdCalled)
	}

	cmd := exec.Command(binPath+"/"+cmdCalled, flag.Args()[1:]...)

	cmd.Stdin = os.Stdin
	cmd.Stdout = os.Stdout
	cmd.Stderr = os.Stderr

	err := cmd.Run()
	exception.FatalOnErr(err)
}
