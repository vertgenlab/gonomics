package main

import (
	"bufio"
	"flag"
	"fmt"
	"github.com/vertgenlab/gonomics/exception"
	"io"
	"io/ioutil"
	"log"
	"os"
	"os/exec"
	"sync"
)

func usage() {
	flag.PrintDefaults()
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

func getGonomicsFuncs() map[string]bool {
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

// from https://gist.github.com/hivefans/ffeaf3964924c943dd7ed83b406bbdea
func streamOutput(output io.ReadCloser, wg *sync.WaitGroup) {
	buf := bufio.NewReader(output)
	num := 1
	for num <= 3 {
		outline, _, _ := buf.ReadLine()
		num += 1
		fmt.Println(string(outline))
	}
	wg.Done()
}

func main() {
	//var list *bool = flag.Bool("list", false, "List cmds available in gonomics.")
	flag.Usage = usage
	log.SetFlags(log.Ldate | log.Ltime | log.Lshortfile)
	flag.Parse()

	if len(flag.Args()) == 0 {
		flag.Usage()
		log.Fatal()
	}

	cmdCalled := flag.Arg(0) // command called by the user (e.g. 'gonomics faFormat')

	binPath, binExists := getBin()
	isGonomicsFunc := getGonomicsFuncs()

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
