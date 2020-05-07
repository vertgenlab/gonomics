package fasta

import(
	"sync"
	//"github.com/vertgenlab/gonomics/common"
	"github.com/vertgenlab/gonomics/dna"
	"github.com/vertgenlab/gonomics/fileio"
	"strings"
)
func GoReadToChan(filename string) <-chan *Fasta {
	output := make(chan *Fasta)
	go ReadToChan(filename, output, nil, true)
	return output
}

func ReadToChan(filename string, output chan<- *Fasta, wg *sync.WaitGroup, closeChannel bool) {
	file := fileio.EasyOpen(filename)
	for fa, done := NextFasta(file); !done; fa, done = NextFasta(file) {
		output <- fa	
	}
	if closeChannel {
		close(output)
	}
	if wg != nil {
    	wg.Done()
    }
}

func FinalContigsToCanu(filterOut string, output string, contigs []string) {
	finalFa := fileio.EasyCreate(output)
	faChan := make(chan *Fasta)
	var wg sync.WaitGroup
	defer finalFa.Close()
	var contamination map[string]string
	if strings.Compare(filterOut, "None") != 0 {
		contamination = toContaminationMap(filterOut)
	} else {
		contamination = make(map[string]string)
	}
	wg.Add(len(contigs))
	for _, file := range contigs {
		go ReadToChan(file, faChan, &wg, false)
	}
	var writer sync.WaitGroup
	writer.Add(1) 
	go FaFilterWrite(finalFa, contamination, faChan, &writer)
	wg.Wait()
	close(faChan)
	writer.Wait()
}

func FaFilterWrite(mergedFile *fileio.EasyWriter, dirty map[string]string, contigs <-chan *Fasta, wg *sync.WaitGroup) {
	for i := range contigs {
		_, key := dirty[i.Name]
		if !key {
			i.Name = strings.Join(strings.Split(i.Name, "_"), " ")
			WriteToFileHandle(mergedFile, i, 50)
		}
	}
	wg.Done()
}

func toContaminationMap(filename string) map[string]string {
	failedFilter := fileio.Read(filename)
	contamination := make(map[string]string)
	var column []string
	for _, blastDb := range failedFilter {
		column = strings.Split(blastDb, "\t")
		//queryName
		//i noticed some of the pacbio adaptors also mapped to other weird bacteria, and wanted to prioritize the pacbio alignment
		//pacbio contracuts are added to the map automatically

		_, ok := contamination[column[0]]
		if !ok {
			contamination[column[0]] = column[1]
		}
		

	}
	return contamination
}

func MultiFileChan(files[]string, output chan<- *Fasta) {
	//output := make(chan *Fasta, 824)
	var wg sync.WaitGroup
	wg.Add(len(files))
	for _, each := range files {
		go ReadToChan(each, output, &wg, false)
	}
	wg.Wait()
	close(output)
}

func NextFasta(reader *fileio.EasyReader) (*Fasta, bool) {
	var fa Fasta
	line, done := fileio.EasyNextLine(reader)
	if done {
		return nil, true
	} else {
		if strings.HasPrefix(line, ">") {
			fa = Fasta{Name: line[1:len(line)], Seq: nextSeq(reader)}
		}
	}
	return &fa, false
}

func nextSeq(reader *fileio.EasyReader) []dna.Base {
	var line string
	var err error
	var nextBytes []byte
	var answer []dna.Base
	for nextBytes, err = reader.Peek(1); len(nextBytes) > 0 && nextBytes[0] != '>' && err == nil; nextBytes, err = reader.Peek(1) {
		line, _ = fileio.EasyNextLine(reader)
		answer = append(answer, dna.StringToBases(line)...)
	}
	return answer
}
/*
func WritingChannelMultiFiles(files []*fileio.EasyWriter, output <-chan *Fasta, wg *sync.WaitGroup) {
	var index int = 0
	for fa := range output {
		WriteToFileHandle(files[index], fa, 50)
		index++
		if index == len(files) {
			index = 0
		}
	}
	wg.Done()
}*/

func WritingChannel(file *fileio.EasyWriter, output <-chan *Fasta, wg *sync.WaitGroup) {
	for fa := range output {
		WriteToFileHandle(file, fa, 50)
	}
	wg.Done()
}


func WriteGroups(filename string, groups [][]*Fasta) {
	file := fileio.EasyCreate(filename)
	defer file.Close()
	lineLength := 50
	for i, _ := range groups {
		for _, records := range groups[i] {
			WriteToFileHandle(file, records, lineLength)
		}
	}
}
