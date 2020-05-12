package chain

import (
	"fmt"
	"sync"
	"testing"
	"github.com/vertgenlab/gonomics/fileio"
)

var filename string = "testdata/gasAcu1ToRabsLiftover.chain"
func TestReader(t *testing.T) {
	chainfile, _ := Read(filename)
	for i := 0; i < len(chainfile); i++ {
		fmt.Print(ChainToString(chainfile[i]))
	}
}

func TestReadToChan(t *testing.T) {
	file := fileio.EasyOpen(filename)
	defer file.Close()
	answer := make(chan *Chain)
	SaveComments(file)
	var wg sync.WaitGroup
	go ReadToChan(file, answer)
	//wg.Add(1)
	go ChanToFile("/dev/stdout", answer, nil, &wg)
}
