// Command Group: "Data Simulation"

// Returns random Pfasta files
package main

import (
	"flag"
	"fmt"
	"github.com/vertgenlab/gonomics/dna/pDna"
	"github.com/vertgenlab/gonomics/fasta/pFasta"
	"log"
	//"math"
	"math/rand"
	//"strconv"
)

func randSeq(length int) pFasta.PFasta {
	var one, two, three, four, total float32
	var answer = pFasta.PFasta{}
	answer.Name = "randSeq"
	answer.Seq = make([]pDna.Float32Base, length)
	for i := 0; i < length; i++ {
		one = rand.Float32()
		two = rand.Float32()
		three = rand.Float32()
		four = rand.Float32()
		total = one + two + three + four
		answer.Seq[i].A = one / total
		answer.Seq[i].C = two / total
		answer.Seq[i].G = three / total
		answer.Seq[i].T = four / total
	}
	return answer
}

func ASeq(length int) pFasta.PFasta {
	var answer = pFasta.PFasta{}
	answer.Name = "ASeq"
	answer.Seq = make([]pDna.Float32Base, length)
	for i := 0; i < length; i++ {
		answer.Seq[i].A = 1
		answer.Seq[i].C = 0
		answer.Seq[i].G = 0
		answer.Seq[i].T = 0
	}
	return answer
}

func ASeq2(length int) pFasta.PFasta {
	var answer = pFasta.PFasta{}
	answer.Name = "ASeq2"
	answer.Seq = make([]pDna.Float32Base, length)
	for i := 0; i < length; i++ {
		answer.Seq[i].A = 1
		answer.Seq[i].C = 0
		answer.Seq[i].G = 0
		answer.Seq[i].T = 0
	}
	return answer
}

func CSeq(length int) pFasta.PFasta {
	var answer = pFasta.PFasta{}
	answer.Name = "CSeq"
	answer.Seq = make([]pDna.Float32Base, length)
	for i := 0; i < length; i++ {
		answer.Seq[i].A = 0
		answer.Seq[i].C = 1
		answer.Seq[i].G = 0
		answer.Seq[i].T = 0
	}
	return answer
}

func GapSeq(length int) pFasta.PFasta {
	var answer = pFasta.PFasta{}
	answer.Name = "GapSeq"
	answer.Seq = make([]pDna.Float32Base, length)
	for i := 0; i < length; i++ {
		answer.Seq[i].A = 0
		answer.Seq[i].C = 0
		answer.Seq[i].G = 0
		answer.Seq[i].T = 0
	}
	return answer
}

// func simulatePfasta(number int, length int) {
func simulatePfasta() {
	var outFile string
	/*
		var answerPfa = make([]pFasta.PFasta, 1) // answerPfa is a [] of pFasta
		var pdna = pDna.Float32Base{
			A: 1,
			C: 0,
			G: 0,
			T: 0,
		}
	*/
	answerPfa := []pFasta.PFasta{randSeq(10000), ASeq(10000), ASeq2(10000), CSeq(10000), GapSeq(10000)} // make answerPfa []pFasta with the randSeq function borrowed from Craig's pFasta Read bug fix
	outFile = "simulatePfasta.pFa"
	pFasta.Write(outFile, answerPfa)
	//var answerPfaRead []pFasta.PFasta

	/*
		for i := 0; i <= number; i++ {
				answerPfa[0].Name = strconv.Itoa(i)
				for j := 0; j <= length; j++ {
					answerPfa[0].Seq = append(answerPfa[0].Seq, pdna)
				}
			outFile = strconv.Itoa(i) + ".pFa"
			pFasta.Write(outFile, answerPfa)
			answerPfaRead = pFasta.Read(outFile)
			for _, v := range answerPfaRead {
				for k := range v.Seq {
					if math.IsNaN(float64(v.Seq[k].A)) || math.IsNaN(float64(v.Seq[k].C)) || math.IsNaN(float64(v.Seq[k].G)) || math.IsNaN(float64(v.Seq[k].T)) || v.Seq[k].A < 0 || v.Seq[k].C < 0 || v.Seq[k].G < 0 || v.Seq[k].T < 0 {
						log.Fatalf("Read pFasta and found NaN or Negative base: %v at position %v, and fasta name is: %v\n", v.Seq[i], i, v.Name)
					}
					if k == 10 || k == 250 || k == 750 {
						fmt.Printf("Read pFasta, Check base: %v at position %v, and fasta name is: %v\n", v.Seq[k], k, v.Name)
					}
				}
				//fmt.Printf("v: %v\n", v)
				fmt.Printf("length of pFasta: %d\n", len(v.Seq))
			}
		}
	*/
}

func usage() {
	fmt.Print(
		"simulatePfasta - Returns random Pfasta files.\n" +
			"Usage:\n" +
			" simulatePfasta number length\n" +
			"options:\n")
	flag.PrintDefaults()
}

func main() {
	//var expectedNumArgs int = 2
	//var number int
	//var length int
	//var err error

	flag.Usage = usage
	log.SetFlags(log.Ldate | log.Ltime | log.Lshortfile)
	flag.Parse()

	/*
		if len(flag.Args()) != expectedNumArgs {
			flag.Usage()
			log.Fatalf("Error: expecting %d arguments, but got %d\n",
				expectedNumArgs, len(flag.Args()))
		}
	*/

	//number, err = strconv.Atoi(flag.Arg(0))
	//length, err = strconv.Atoi(flag.Arg(1))
	/*
		if err != nil {
			panic(err)
		}
	*/

	//simulatePfasta(number, length)
	simulatePfasta()
}
