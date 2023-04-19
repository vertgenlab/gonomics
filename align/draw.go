package align

import (
	"fmt"
	"image"
	"image/color"
	"image/draw"
	"sort"
	"strings"

	"github.com/vertgenlab/gonomics/dna"
	"github.com/vertgenlab/gonomics/fasta"
	"github.com/vertgenlab/gonomics/sketch"
)

type keyValue struct {
	Key   string
	Value int
}

func incrementOrAdd(list []keyValue, needle string) []keyValue {
	for i := range list {
		if list[i].Key == needle {
			list[i].Value++
			return list
		}
	}
	return append(list, keyValue{Key: needle, Value: 1})
}

func determineChunkColors(aln []fasta.Fasta, chunkSize int, palette color.Palette) (map[string]color.Color, error) {
	answer := make(map[string]color.Color, 0)
	list := make([]keyValue, 0)

	for i := range aln {
		if len(aln[i].Seq)%chunkSize != 0 {
			return nil, fmt.Errorf("The %s sequence has a length of %d, which is not divisible by a chunkSize of %d\n", aln[i].Name, len(aln[i].Seq), chunkSize)
		}
		for chunkStart := 0; chunkStart < len(aln[i].Seq); chunkStart += chunkSize {
			chunkText := dna.BasesToString(aln[i].Seq[chunkStart:(chunkStart + chunkSize)])
			gapCount := strings.Count(chunkText, "-")
			switch gapCount {
			case chunkSize: /* ignore if all gaps */
			case 0:
				list = incrementOrAdd(list, chunkText)
			default:
				return nil, fmt.Errorf("Error: %s should be either all gaps or no gaps\n", chunkText)
			}
		}
	}
	fmt.Printf("Number of chunks recorded: %d\n", len(list))
	sort.Slice(list, func(i, j int) bool { return list[i].Value > list[j].Value })

	for i := 0; i < len(list) && i < len(palette); i++ {
		answer[list[i].Key] = palette[i]
	}

	return answer, nil
}

func drawNames(aln []fasta.Fasta, vSpacing int) *image.RGBA {
	imageWidth := 200
	imageHeight := len(aln) * vSpacing
	img := image.NewRGBA(image.Rect(0, 0, imageWidth, imageHeight))
	sketch.FilledRectangle(img, 0, 0, imageWidth, imageHeight, color.White) /* make everything white to start */
	for i := range aln {
		sketch.Text(img, aln[i].Name, 0, (i+1)*vSpacing)
	}
	return img
}

func DrawAlignedChunks(aln []fasta.Fasta, chunkSize int, chunkPixelWidth int, chunkPixelHeight int) (*image.RGBA, error) {
	imgChunks, err := drawChunks(aln, chunkSize, chunkPixelWidth, chunkPixelHeight)
	if err != nil {
		return nil, err
	}
	imgNames := drawNames(aln, chunkPixelHeight)
	imageHeight := imgChunks.Bounds().Max.Y - imgChunks.Bounds().Min.Y
	imageWidth := imgChunks.Bounds().Max.X - imgChunks.Bounds().Min.X + 10 + imgNames.Bounds().Max.X - imgNames.Bounds().Min.X
	img := image.NewRGBA(image.Rect(0, 0, imageWidth, imageHeight))
	sketch.FilledRectangle(img, 0, 0, imageWidth, imageHeight, color.White) /* make everything white to start */
	draw.Draw(img, imgChunks.Bounds(), imgChunks, image.Point{X: 0, Y: 0}, draw.Src)
	rect := imgNames.Bounds()
	rect.Min.X += imgChunks.Bounds().Max.X - imgChunks.Bounds().Min.X + 10
	rect.Max.X += imgChunks.Bounds().Max.X - imgChunks.Bounds().Min.X + 10
	draw.Draw(img, rect, imgNames, image.Point{X: 0, Y: 0}, draw.Src)
	return img, nil
}

func drawChunks(aln []fasta.Fasta, chunkSize int, chunkPixelWidth int, chunkPixelHeight int) (*image.RGBA, error) {
	colorMap, err := determineChunkColors(aln, chunkSize, sketch.TrubetskoyPalette[:19])
	if err != nil {
		return nil, err
	}
	allGaps := strings.Repeat("-", chunkSize)
	colorMap[allGaps] = color.Black

	alnLength := len(aln[0].Seq)
	numSeq := len(aln)
	imageWidth := alnLength / chunkSize * chunkPixelWidth
	imageHeight := chunkPixelHeight * numSeq
	img := image.NewRGBA(image.Rect(0, 0, imageWidth, imageHeight))
	sketch.FilledRectangle(img, 0, 0, imageWidth, imageHeight, color.White) /* make everything white to start */
	for i := range aln {
		for chunkStart := 0; chunkStart < len(aln[i].Seq); chunkStart += chunkSize {
			chunkText := dna.BasesToString(aln[i].Seq[chunkStart:(chunkStart + chunkSize)])
			chunkColor, found := colorMap[chunkText]
			if !found {
				chunkColor = sketch.TrubetskoyPalette[19] /* gray */
			}
			xStart := chunkStart / chunkSize * chunkPixelWidth
			xEnd := xStart + chunkPixelWidth
			yStart := i * chunkPixelHeight
			yEnd := yStart + chunkPixelHeight
			sketch.FilledRectangle(img, xStart, yStart, xEnd, yEnd, chunkColor)
		}
	}
	return img, nil
}
