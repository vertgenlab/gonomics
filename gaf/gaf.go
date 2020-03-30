package gaf

import (
//"github.com/vertgenlab/gonomics/common"
)

type Gaf struct {
	QName    string
	QLen     int64
	QStart   int64
	QEnd     int64
	QStrand  bool
	Path     *Path
	PLen     int64
	PStart   int64
	PEnd     int64
	MatchNum uint32
	BlkLen   uint32
	MapQ     int32
	Notes    string
}

//TODO: Figure out how to define stable path...
type Path struct {
	Name       string
	StablePath bool
	Start      uint32
	End        uint32
	Next       *Path
}

//TODO:
/*
func Read(filename string) []Gaf {
	file := fileio.EasyOpen(filename)
	var line string
	var doneReading bool = false
	var text []string
	defer file.Close()
	for line, doneReading = fileio.EasyNextRealLine(file); !doneReading; line, doneReading = fileio.EasyNextRealLine(file) {


	}
}

func gafLine(line string) *Gaf {
	var curr Gaf
	text = strings.SplitN(line, "\t", 13)
	if len(text) < 12 {
		log.Fatal(fmt.Errorf("Was expecting atleast 12 columns per line, but this line did not:%s\n", line))
	}
	curr.QName = text[0]
	curr.QLen = common.StringToInt64(text[1])
	curr.QStart = common.StringToInt64(text[2])
	curr.QEnd = common.StringToInt64(text[3])
	curr.QStrand = text[4]
	curr.Path = common.StringToInt64(text[5]
	curr.PLen = common.StringToInt64(text[6])
	curr.PStart = common.StringToInt64(text[7])
	curr.PEnd = common.StringToInt64(text[8])
	curr.MatchNum = common.StringToUint32(text[9])
	curr.BlkLen = common.StringToUint32(text[10])
	curr.MapQ = common.StringToInt64(text[11])
	if len(text) == 13 {
		curr.Notes = text[12]
	}
}

func pathToString(tab string) *Path {

}*/
