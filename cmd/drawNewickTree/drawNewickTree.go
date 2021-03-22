package main

import (
	"flag"
	"fmt"
	"github.com/vertgenlab/gonomics/tree"
	"image"
	"image/png"
	"log"
	"os"
)

//drawNewickTree takes in a newick file and returns a png of the newick tree, where the image's width and height are as specified, or default to 1500
func drawNewickTree(newickFile string, pngFile string, imgWidth int, imgHeight int) {
	var nt *tree.Tree
	var err error
	var img *image.RGBA
	var imgOutFile *os.File

	nt, err = tree.ReadNewick(newickFile)
	if err != nil {
		log.Panic(err)
	}
	img, err = tree.Draw(nt, imgWidth, imgHeight)
	if err != nil {
		log.Panicf("Error in Draw")
	}
	imgOutFile, err = os.Create(pngFile)
	if err != nil {
		log.Panicf("Error in Create")
	}
	defer imgOutFile.Close()
	err = png.Encode(imgOutFile, img)
	if err != nil {
		log.Panicf("Error in png.Encode")
	}
}

func usage() {
	fmt.Print(
		"drawNewickTree takes in a newick format text file and outputs a png for tree visualization. If image width and height aren't specified they will default to 1500\n" +
			"Usage:\n" +
			"drawNewicktree [-option int] <input.txt> <output.png>\n" +
			"options:\n")
	flag.PrintDefaults()
}

func main() {
	var expectedNumArgs = 2
	var imgWidth *int = flag.Int("imgWidth", 1500, "Specifies the width of the ouput.png.")
	var imgHeight *int = flag.Int("imgHeight", 1500, "Specifies the height of the ouput.png.")

	flag.Usage = usage
	log.SetFlags(log.Ldate | log.Ltime | log.Lshortfile)
	flag.Parse()

	if len(flag.Args()) != expectedNumArgs {
		flag.Usage()
		log.Fatalf("Error: expecting %d arguments, but got %d\n",
			expectedNumArgs, len(flag.Args()))
	}

	infile := flag.Arg(0)
	imgOutFile := flag.Arg(1)

	drawNewickTree(infile, imgOutFile, *imgWidth, *imgHeight)
}
