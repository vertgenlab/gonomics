package vcf

import (
	"bufio"
	"github.com/vertgenlab/gonomics/common"
	"io"
	"os"
	"strings"
	"sync"
	"log"
)

func ReadToChan(filename string) chan *Vcf {
	var answer = make(chan *Vcf)
	var wg sync.WaitGroup
	wg.Add(1)

	go func() {
		var curr *Vcf
		var line string
		file, err := os.Open(filename)
		if err != nil {
			log.Fatalln("error opening input file")
		}
		reader := bufio.NewReader(file)
		if err != nil {
			log.Fatalln("error reading input file")
		}
		var err2 error
		//var rline []byte
		for ; err2 != io.EOF; line, err2 = reader.ReadString('\n') {
			//line = string(rline[:])
			data := strings.Split(line, "\t")
			//fmt.Println("there is data here")
			switch {
			case strings.HasPrefix(line, "#"):
				//don't do anything
			case len(data) == 1:
				//these lines are sequences, and we are not recording them
				//fmt.Println("found sequences")
			case len(line) == 0:
				//blank line

			case len(data) == 10:
				curr = &Vcf{Chr: data[0], Pos: common.StringToInt64(data[1]), Id: data[2], Ref: data[3], Alt: data[4], Qual: common.StringToFloat64(data[5]), Filter: data[6], Info: data[7], Format: data[8], Notes: data[9]}
				answer <- curr
			default:
				//fmt.Println("unexpected line")
			}
		}
		wg.Done()
	}()

	go func() {
		wg.Wait()
		close(answer)
	}()

	return answer
}
