package sketch

import (
	"image"
	"image/color"
	"testing"
)

// TestHLine tests the HLine function by drawing a horizontal line and checking if the expected color is set at each pixel.
func TestHLine(t *testing.T) {
	img := image.NewRGBA(image.Rect(0, 0, 10, 10))
	col := color.RGBA{255, 0, 0, 255}

	HLine(img, 2, 7, 5, col)

	for x := 2; x < 7; x++ {
		// This test ensures that the HLine function correctly sets the color for a horizontal line segment on the image.
		if got := img.At(x, 5); got != col {
			t.Errorf("Error: Expected color %v at position (%d, %d), got %v", col, x, 5, got)
		}
	}
}

// TestVLine tests the VLine function by drawing a vertical line and checking if the expected color is set at each pixel.
func TestVLine(t *testing.T) {
	img := image.NewRGBA(image.Rect(0, 0, 10, 10))
	col := color.RGBA{0, 255, 0, 255}

	VLine(img, 5, 2, 7, col)
	// This test ensures that the VLine function correctly sets the color for a vertical line segment on the image.
	for y := 2; y < 7; y++ {
		if got := img.At(5, y); got != col {
			t.Errorf("Error: Expected color %v at position (%d, %d), got %v", col, 5, y, got)
		}
	}
}

// TestRectangle tests the Rectangle function by drawing a rectangle and checking if the expected color is set at each pixel.
func TestRectangle(t *testing.T) {
	img := image.NewRGBA(image.Rect(0, 0, 10, 10))
	col := color.RGBA{0, 0, 255, 255}

	Rectangle(img, 2, 2, 7, 7, col)
	// This test ensures that the Rectangle function correctly sets the color for all the pixels within the rectangle on the image.
	for x := 2; x < 7; x++ {
		if got := img.At(x, 2); got != col {
			t.Errorf("Error: Expected color %v at position (%d, %d), got %v", col, x, 2, got)
		}
		if got := img.At(x, 7); got != col {
			t.Errorf("Error: Expected color %v at position (%d, %d), got %v", col, x, 7, got)
		}
	}

	for y := 2; y < 7; y++ {
		if got := img.At(2, y); got != col {
			t.Errorf("Error: Expected color %v at position (%d, %d), got %v", col, 2, y, got)
		}
		if got := img.At(7, y); got != col {
			t.Errorf("Error: Expected color %v at position (%d, %d), got %v", col, 7, y, got)
		}
	}
}

// TestFilledRectangle tests the FilledRectangle function by drawing a filled rectangle and checking if the expected color is set at each pixel.
func TestFilledRectangle(t *testing.T) {
	img := image.NewRGBA(image.Rect(0, 0, 10, 10))
	col := color.RGBA{255, 0, 255, 255}

	FilledRectangle(img, 2, 2, 7, 7, col)
	// This test ensures that the FilledRectangle function correctly sets the color for all the pixels within the filled rectangle on the image.
	for x := 2; x < 7; x++ {
		for y := 2; y < 7; y++ {
			if got := img.At(x, y); got != col {
				t.Errorf("Error: Expected color %v at position (%d, %d), got %v", col, x, y, got)
			}
		}
	}
}
