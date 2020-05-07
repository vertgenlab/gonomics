package main

import (
	"fmt"
	"github.com/vertgenlab/gonomics/common"
	"github.com/vertgenlab/gonomics/fasta"
	"github.com/vertgenlab/gonomics/fileio"
	"log"
	"strings"
	"sync"
)

func taxonIdsToMap(filename string) map[int]bool {
	taxChan := make(chan string)
	go toChan(filename, taxChan)
	db := make(map[int]bool)
	var curr int
	for line := range taxChan {
		curr = common.StringToInt(line)
		_, ok := db[curr]
		if !ok {
			db[curr] = true
		}
	}
	return db
}

func filterMatchingTax(db map[int]bool, line string) bool {
	columns := strings.Split(line, "\t")
	id := common.StringToInt(columns[4])
	return db[id]
}

func FinalContigsToCanu(filterOut string, output string, contigs []string) {
	finalFa := fileio.EasyCreate(output)

	faChan := make(chan *fasta.Fasta)
	var wg sync.WaitGroup
	defer finalFa.Close()
	var dbFilter map[string]bool
	//input none if you do not want anything filtered, and just want fasta merged
	if strings.Compare(filterOut, "None") == 0 {
		dbFilter = make(map[string]bool)
	} else {
		dbFilter = blastMap(filterOut)
	}
	wg.Add(len(contigs))
	for _, file := range contigs {
		go fasta.ReadToChan(file, faChan, &wg, false)
	}
	var writer sync.WaitGroup
	writer.Add(1)
	go FaFilterWrite(finalFa, dbFilter, faChan, &writer)
	wg.Wait()
	close(faChan)
	writer.Wait()
}

func FaFilterWrite(mergedFile *fileio.EasyWriter, db map[string]bool, contigs <-chan *fasta.Fasta, wg *sync.WaitGroup) {
	for eachRead := range contigs {
		_, key := db[eachRead.Name]
		if !key {
			eachRead.Name = strings.Join(strings.Split(eachRead.Name, "_"), " ")
			fasta.WriteToFileHandle(mergedFile, eachRead, 50)
		}
	}
	wg.Done()
}

//takes a 2 couln list (read name, target blastHit) of read names and makes a map of reads to filter to filter sequence reads
func blastMap(filename string) map[string]bool {
	failedFilter := fileio.Read(filename)
	db := make(map[string]bool)
	var column []string
	for _, read := range failedFilter {
		column = strings.Split(read, "\t")
		_, ok := db[column[0]]
		if !ok {
			db[column[0]] = true
		}
	}
	return db
}

func WrapAllDbFiles(output string, taxons string, files []string) {
	final := fileio.EasyCreate(output)
	defer final.Close()
	progress := fileio.EasyCreate("progress.txt")
	defer progress.Close()
	txDb := taxonIdsToMap(taxons)
	for _, each := range files {
		writeSummary(progress, fmt.Sprintf("%s", each))
		readText(each, txDb, final)
	}
}

func toChan(filename string, output chan<- string) {
	file := fileio.EasyOpen(filename)
	defer file.Close()
	for line, doneReading := fileio.EasyNextLine(file); !doneReading; line, doneReading = fileio.EasyNextLine(file) {
		output <- line
	}
	close(output)
}

//queryName, quertStart, queryEnd, targetName, targetStart, targetEnd, (6)score, bitScore, (8)alignLength, eValue, percentAlign, numOfMatches, mismatch, querySeq, targetSeq
func readText(filename string, txDb map[int]bool, output *fileio.EasyWriter) {
	blast := make(chan string)
	go toChan(filename, blast)
	for text := range blast {
		if !filterMatchingTax(txDb, text) {
			writeSummary(output, text)
			log.Printf("%s\n", text)
		}
	}
}

func writeSummary(file *fileio.EasyWriter, failedFilter string) {
	var err error
	_, err = fmt.Fprintf(file, "%s\n", failedFilter)
	common.ExitIfError(err)
}

/*
func sticklebackHit(subject string) bool {
	if strings.Contains(subject, "Gasterosteus") || strings.Contains(subject, "aculeatus") || strings.Contains(subject, "Pungitius") {
		return true
	}
	if strings.Contains(subject, "Apeltes") || strings.Contains(subject, "quadracus") {
		return true
	}
	if strings.Contains(subject, "Culaea") || strings.Contains(subject, "inconstans") {
		return true
	}
	return false
}

func contaminationDb(subject string) bool {
	db := []string{"Synthetic_construct_PacBio_DNA_internal_control_sequence",
		"Synthetic_construct_PacBio_unrolled_DNA_internal_control_sequence",
	}
	for _, hit := range db {
		if strings.Contains(subject, hit) {
			return true
		}
	}
	return false
}

func rnaFilter(subject string) bool {
	if strings.Contains(subject, "mRNA") {
		return false
	}
	if strings.Contains(subject, "ribosomal_RNA") || strings.Contains(subject, "rRNA") {
		return false
	}
	if strings.Contains(subject, "partial_cds") {
		return false
	}
	return true
}

func alignmentLengthFilter(query string, alignLen string) bool {
	qLen := common.StringToFloat64(strings.Split(query, ",")[1])
	if common.StringToFloat64(alignLen)/qLen < 0.5 {
		return false
	} else {
		return true
	}
}*/
