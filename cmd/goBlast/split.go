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

func FaToMatrix(fa []*fasta.Fasta, numSplit int) [][]*fasta.Fasta {
	var i, j int
	var faidx int
	groups := len(fa) / numSplit
	leftover := len(fa) % numSplit

	matrix := make([][]*fasta.Fasta, groups)
	log.Printf("Total number of contigs/scaffolds: %d\n", len(fa))
	log.Printf("Divide into %d groups containing %d each and %d leftover\n", numSplit, groups, leftover)
	for column := range matrix {
		matrix[column] = make([]*fasta.Fasta, numSplit)
	}
	for i = 0; i < groups; i++ {
		for j = 0; j < numSplit; j++ {
			matrix[i][j] = fa[faidx]
			faidx++
		}
	}

	for i = 0; i < leftover; i++ {
		matrix[i] = append(matrix[i], fa[faidx])
		faidx++
	}

	if faidx != len(fa) {
		log.Fatalf("Error: Number of elements in fasta matrix do not match input fasta")
	}
	return matrix
}

/*
func TransposeFaMatrix(matrix [][]*Fasta) [][]*Fasta {
	xl := len(matrix[0])
	yl := len(matrix)
	leftover := len(matrix[0]) - len(matrix[len(matrix)-1])
	result := make([][]*Fasta, xl)
	for i := range result {
		result[i] = make([]*Fasta, yl)
		for k := 0; k < len(result[i]); k++ {
			result[i][k] = &Fasta{}
		}
	}
	for i := 0; i < xl-leftover; i++ {
		for j := 0; j < yl; j++ {

			result[i][j] = matrix[j][i]
		}
	}
	return result
}*/

func TransposeFaMatrix(matrix [][]*fasta.Fasta) [][]*fasta.Fasta {
	//x
	rows := len(matrix[0]) - 1 //handle possible leftovers from uneven matrices
	//y
	columns := len(matrix)
	transposed := make([][]*fasta.Fasta, rows)
	var i, j, k int
	for i = 0; i < len(transposed); i++ {
		transposed[i] = make([]*fasta.Fasta, columns)
		for j = 0; j < len(transposed[i]); j++ {
			transposed[i][j] = &fasta.Fasta{}
		}
	}
	for i = 0; i < rows; i++ {
		for j = 0; j < columns; j++ {
			transposed[i][j] = matrix[j][i]
		}
	}
	extra := len(matrix[0]) - len(matrix[len(matrix)-1])
	if extra > 0 {
		var leftovers []*fasta.Fasta

		for k = 0; k < columns; k++ {
			if len(matrix[k]) > rows {
				leftovers = append(leftovers, matrix[k][rows])
			}
		}
		transposed = append(transposed, leftovers)
	}
	return transposed
}

func EfficientTranspose(filename string, prefix string, groups int) {
	faChan := make(chan *fasta.Fasta)
	go fasta.ReadToChan(filename, faChan, nil, true)
	splitFiles := make([]*fileio.EasyWriter, groups)
	var err error
	for i := 0; i < len(splitFiles); i++ {
		splitFiles[i] = fileio.EasyCreate(fmt.Sprintf("%s_%d.fa", prefix, i))
		defer splitFiles[i].Close()
		common.ExitIfError(err)
	}
	var wg sync.WaitGroup
	wg.Add(1)
	go WritingChannelMultiFiles(splitFiles, faChan, &wg)
	wg.Wait()
}

func WritingChannelMultiFiles(files []*fileio.EasyWriter, output <-chan *fasta.Fasta, wg *sync.WaitGroup) {
	var index int = 0
	for fa := range output {
		fa.Name = strings.Join(strings.Split(fa.Name, " "), "_")
		fasta.WriteToFileHandle(files[index], fa, 50)
		index++
		if index == len(files) {
			index = 0
		}
	}
	wg.Done()
}
