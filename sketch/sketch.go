// Package sketch provides utilities and color palettes for plotting graphs and images
package sketch

import (
	"image"
	"image/color"

	"golang.org/x/image/font"
	"golang.org/x/image/font/basicfont"
	"golang.org/x/image/math/fixed"
)

// HLine draws a horizontal line on the given image between the xStart and xEnd coordinates at the y coordinate with the specified color.
func HLine(img *image.RGBA, xStart, xEnd, y int, col color.Color) {
	for x := xStart; x < xEnd; x++ {
		img.Set(x, y, col)
	}
}

// VLine draws a vertical line on the given image between the yStart and yEnd coordinates at the x coordinate with the specified color.
func VLine(img *image.RGBA, x, yStart, yEnd int, col color.Color) {
	for y := yStart; y < yEnd; y++ {
		img.Set(x, y, col)
	}
}

// Rectangle draws a rectangle on the given image with the specified coordinates and color.
func Rectangle(img *image.RGBA, xOne, yOne, xTwo, yTwo int, col color.Color) {
	HLine(img, xOne, xTwo, yOne, col)
	HLine(img, xOne, xTwo, yTwo, col)
	VLine(img, xOne, yOne, yTwo, col)
	VLine(img, xTwo, yOne, yTwo, col)
}

// FilledRectangle draws a filled rectangle on the given image with the specified coordinates and color.
func FilledRectangle(img *image.RGBA, xOne, yOne, xTwo, yTwo int, col color.Color) {
	for x := xOne; x < xTwo; x++ {
		for y := yOne; y < yTwo; y++ {
			img.Set(x, y, col)
		}
	}
}

// Text writes the specified label on the given image at the xStart and yStart coordinates.
func Text(img *image.RGBA, label string, xStart, yStart int) {
	point := fixed.Point26_6{fixed.Int26_6(xStart * 64), fixed.Int26_6(yStart * 64)}

	d := &font.Drawer{
		Dst:  img,
		Src:  image.NewUniform(color.Black),
		Face: basicfont.Face7x13,
		// Face: inconsolata.Regular8x16,
		Dot: point,
	}
	d.DrawString(label)
}
