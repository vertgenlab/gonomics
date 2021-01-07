//template copied from mafFilter.go
package main

import (
	"flag"
	"fmt"
	"github.com/vertgenlab/gonomics/maf"
	"log"
)

//template copied from Read function in maf.go
func Read(filename string) []*Maf {
	var answer []*Maf
	var line, prevLine string
	var doneReading bool = false
	var words []string
	var curr *Maf
	var currSpecies *MafSpecies

	file := fileio.EasyOpen(filename)
	defer file.Close()

	for line, doneReading = fileio.EasyNextRealLine(file); !doneReading; line, doneReading = fileio.EasyNextRealLine(file) {
		if strings.HasPrefix(line, "a") {
			if curr != nil {
				log.Fatalf("Error: no blank line before another 'a' line at line: %s\n", line)
			}
			curr = parseMafALine(line)
		} else if strings.HasPrefix(line, "s") || strings.HasPrefix(line, "i") || strings.HasPrefix(line, "e") {
			if curr == nil {
				log.Fatalf("Error: did not find an 'a' line before this, 'sie' line: %s\n", line)
			}
			words = strings.Fields(line)
			currSpecies = FindSpeciesExactMatch(curr, words[1])
			if currSpecies == nil {
				currSpecies = &MafSpecies{Src: words[1]}
				curr.Species = append(curr.Species, currSpecies)
			}
			if strings.HasPrefix(line, "s") {
				if currSpecies.SLine != nil {
					log.Fatalf("Error: this 's' line looks like a duplicate: %s\n", line)
				}
				currSpecies.SLine = parseMafSLine(line)
			} else if strings.HasPrefix(line, "i") {
				if currSpecies.ILine != nil {
					log.Fatalf("Error: this 'i' line looks like a duplicate: %s\n", line)
				}
				currSpecies.ILine = parseMafILine(line)
			} else if strings.HasPrefix(line, "e") {
				if currSpecies.ELine != nil {
					log.Fatalf("Error: this 'e' line looks like a duplicate: %s\n", line)
				}
				currSpecies.ELine = parseMafELine(line)
			} else {
				log.Fatalf("Error: trouble parsing maf line: %s\n", line)
			}
		} else if line == "" { //blank line at end of maf block
			answer = append(answer, curr)
			curr = nil
		} else {
			log.Fatalf("Unexpected format in maf file on line: %s\n", line)
		}
		prevLine = line
	}
	if prevLine != "" {
		log.Fatalf("Error: maf should have a blank line as the last non-comment line, but found this at end: %s\n", prevLine)
	}
	return answer
}

func mafFilter(inFile string, outFile string, threshold float64) {
	mafRecords := maf.Read(inFile)
	var outMaf []*maf.Maf

	for i, _ := range mafRecords {
		if mafRecords[i].Score >= threshold {
			outMaf = append(outMaf, mafRecords[i])
		}
	}

	maf.Write(outFile, outMaf)
}

func usage() {
	fmt.Print(
		"mafFilter - Filter a maf file to remove entries below a score threshold\n" +
			"Usage:\n" +
			" mafFilter mafFile oufMaf\n" +
			"options:\n")
	flag.PrintDefaults()
}

func main() {
	//raven started work here
	//refer to maf/maf.go
	//read file line by line, use maf.Read
	if strings.HasPrefix(line,"e") { //line string
		MafELine=parseMafELine(line)
		if MafELine.Status=='C' && MafELine.Src=="panTro6" {//eC lines, which should belong to Chimp
			//add to collection
		}
	}

	var expectedNumArgs int = 2

	flag.Usage = usage
	var threshold *float64 = flag.Float64("threshold", 0, "Specifies the threshold value")
	log.SetFlags(log.Ldate | log.Ltime | log.Lshortfile)
	flag.Parse()

	if len(flag.Args()) != expectedNumArgs {
		flag.Usage()
		log.Fatalf("Error: expecting %d arguments, but got %d\n", expectedNumArgs, len(flag.Args()))
	}

	mafFile := flag.Arg(0)
	outMaf := flag.Arg(1)

	mafFilter(mafFile, outMaf, *threshold)
}
