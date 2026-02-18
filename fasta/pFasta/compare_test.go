package pFasta

import (
	"github.com/vertgenlab/gonomics/dna/pDna"
    "github.com/vertgenlab/gonomics/dna"
    "github.com/vertgenlab/gonomics/fasta"
    "github.com/vertgenlab/gonomics/wig"
	"testing"
)

// DistTrackTests tests a valid input to DistTrack
var DistTrackTests = []struct {
	Input1     PFasta
    Input2     PFasta
    OutName     string
    DefaultValue    float64
    Expected    map[string]wig.Wig
    Precision   float64
}{
	{Input1: PFasta{Name: "chr1", // first test uses two identical pFas, ensuring the output wig is all zeroes
            Seq: []pDna.Float32Base{
                {A: 0.2,  C: 0.3,  G: 0.4,  T: 0.1},
                {A: 0.25, C: 0.25, G: 0.25, T: 0.25},
                {A: 0.2,  C: 0.3,  G: 0.3,  T: 0.2},
                {A: 0.1,  C: 0.2,  G: 0.3,  T: 0.4},
                {A: 0.6,  C: 0.2,  G: 0.1,  T: 0.1},
            },
		},
		Input2: PFasta{Name: "chr1",
            Seq: []pDna.Float32Base{
                {A: 0.2,  C: 0.3,  G: 0.4,  T: 0.1},
                {A: 0.25, C: 0.25, G: 0.25, T: 0.25},
                {A: 0.2,  C: 0.3,  G: 0.3,  T: 0.2},
                {A: 0.1,  C: 0.2,  G: 0.3,  T: 0.4},
                {A: 0.6,  C: 0.2,  G: 0.1,  T: 0.1},
            },
		},
		OutName:      "chr1.dist",
		DefaultValue: 0.0,
		Expected: map[string]wig.Wig{
			"chr1": wig.Wig{StepType:     "fixedStep",
			Chrom:        "chr1.dist",
			Start:        1,
			Step:         1,
			Span:         1,
			DefaultValue: 0.0,
			Values: []float64{0, 0, 0, 0, 0},},
		},
        Precision: 0.001,
	},
    {Input1: PFasta{Name: "chr1",
            Seq: []pDna.Float32Base{
                {A: 0.2,  C: 0.3,  G: 0.4,  T: 0.1},
                {A: 0.25, C: 0.25, G: 0.25, T: 0.25},
                {A: 0.2,  C: 0.3,  G: 0.3,  T: 0.2},
                {A: 0.1,  C: 0.2,  G: 0.3,  T: 0.4},
                {A: 0.6,  C: 0.2,  G: 0.1,  T: 0.1},
            },
		},
		Input2: PFasta{Name: "chr1",
            Seq: []pDna.Float32Base{
                {A: 0.2,  C: 0.4,  G: 0.3,  T: 0.1},
                {A: 0.15, C: 0.35, G: 0.45, T: 0.05},
                {A: 0.3,  C: 0.2,  G: 0.4,  T: 0.1},
                {A: 0.0,  C: 0.0,  G: 0.5,  T: 0.5},
                {A: 0.0,  C: 0.0,  G: 0.0,  T: 1.0},
            },
		},
		OutName:      "chr1.dist",
		DefaultValue: 0.0,
		Expected: map[string]wig.Wig{
			"chr1": wig.Wig{StepType:     "fixedStep",
			Chrom:        "chr1.dist",
			Start:        1,
			Step:         1,
			Span:         1,
			DefaultValue: 0.0,
			Values: []float64{0.14142136, 0.31622777, 0.2, 0.31622777, 1.10453610},},
		},
        Precision: 0.001,
	},
}


func TestDistTrack(t *testing.T) {
    var res map[string]wig.Wig
	for _, v := range DistTrackTests {
		res = map[string]wig.Wig{"chr1": DistTrack(v.Input1, v.Input2, v.OutName, v.DefaultValue)}
		if !wig.AllEqual(res, v.Expected, v.Precision) {
			t.Errorf("Error: in pFasta. DistTrack valid input test was not as expected.")
		}
        // Write("testdata/test_distTrack_input_A_2.pfa", []PFasta{v.Input1})
        // Write("testdata/test_distTrack_input_B_2.pfa", []PFasta{v.Input2})
	}
}



// DistTrackFastaTests tests a valid input to DistTrackFasta
var DistTrackFastaTests = []struct {
	Input1     PFasta
    Input2     fasta.Fasta
    OutName     string
    DefaultValue    float64
    Expected    map[string]wig.Wig
    Precision   float64
}{
	{Input1: PFasta{Name: "chr1",
            Seq: []pDna.Float32Base{
                {A: 0.2,  C: 0.3,  G: 0.4,  T: 0.1},
                {A: 0.25, C: 0.25, G: 0.25, T: 0.25},
                {A: 0.2,  C: 0.3,  G: 0.3,  T: 0.2},
                {A: 0.1,  C: 0.2,  G: 0.3,  T: 0.4},
                {A: 0.6,  C: 0.2,  G: 0.1,  T: 0.1},
            },
		},
		Input2: fasta.Fasta{Name: "chr1", Seq: dna.StringToBases("TGTAG")},
		OutName:      "chr1.dist",
		DefaultValue: 0.0,
		Expected: map[string]wig.Wig{
			"chr1.dist": wig.Wig{StepType:     "fixedStep",
			Chrom:        "chr1.dist",
			Start:        1,
			Step:         1,
			Span:         1,
			DefaultValue: 0.0,
			Values: []float64{1.048808833,  0.86602540378, 0.927361868188, 1.04880883396, 1.1045360959176},},
		},
        Precision: 0.001,
	},
    {Input1: PFasta{Name: "chr1_map",
            Seq: []pDna.Float32Base{
                {A: 0.2,  C: 0.4,  G: 0.3,  T: 0.1}, // C
                {A: 0.15, C: 0.35, G: 0.45, T: 0.05}, // G
                {A: 0.3,  C: 0.2,  G: 0.4,  T: 0.1}, // G
                {A: 0.0,  C: 0.0,  G: 0.5,  T: 0.5}, // G
                {A: 0.0,  C: 0.0,  G: 0.0,  T: 1.0}, // A
            },
		},
        Input2: fasta.Fasta{Name: "chr1_map", Seq: dna.StringToBases("CGGGA")},
		OutName:      "chr1_map.dist",
		DefaultValue: 0.0,
		Expected: map[string]wig.Wig{
			"chr1_map.dist": wig.Wig{StepType:     "fixedStep",
			Chrom:        "chr1_map.dist",
			Start:        1,
			Step:         1,
			Span:         1,
			DefaultValue: 0.0,
			Values: []float64{0.707107, 0.67082, 0.707107, 0.707107, 1.41421356},},
		},
        Precision: 0.001,
	},
}


func TestDistTrackFasta(t *testing.T) {
    var res map[string]wig.Wig
	for _, v := range DistTrackFastaTests {
		res = map[string]wig.Wig{v.OutName: DistTrackFasta(v.Input1, v.Input2, v.OutName, v.DefaultValue)}
		if !wig.AllEqual(res, v.Expected, v.Precision) {
            t.Errorf("Error: in pFasta. DistTrackFasta valid input test was not as expected.")
		}
        // Write("testdata/test_distTrackFasta_input_A_1.pfa", []PFasta{v.Input1}) // these Write commands are used to generate input files for the cmd/pFaTools distTrack subcommand testing
        // fasta.Write("testdata/test_distTrackFasta_input_B_1.fa", []fasta.Fasta{v.Input2})	
	}
}



// // DistTrackMultiTests tests a valid input to DistTrack
// var DistTrackMultiTests = []struct {
// 	Input1     []PFasta
//     Name1  string
//     Input2     []PFasta
//     Name2   string
//     OutName     string
//     DefaultValue    float64
//     Expected    map[string]wig.Wig
//     Precision   float64
// }{
// 	{
//         Input1: []PFasta{
//             PFasta{Name: "chr1",
//                 Seq: []pDna.Float32Base{
//                     {A: 0.2,  C: 0.3,  G: 0.4,  T: 0.1},
//                     {A: 0.25, C: 0.25, G: 0.25, T: 0.25},
//                     {A: 0.2,  C: 0.3,  G: 0.3,  T: 0.2},
//                     {A: 0.1,  C: 0.2,  G: 0.3,  T: 0.4},
//                     {A: 0.6,  C: 0.2,  G: 0.1,  T: 0.1},},
//                 },
//             PFasta{Name: "chr2",
//                 Seq: []pDna.Float32Base{
//                     {A: 0.5,  C: 0.5,  G: 0.0,  T: 0.0},
//                     {A: 0.5,  C: 0.5,  G: 0.0,  T: 0.0},
//                     {A: 0.5,  C: 0.5,  G: 0.0,  T: 0.0},
//                     {A: 0.5,  C: 0.5,  G: 0.0,  T: 0.0},
//                     {A: 0.5,  C: 0.5,  G: 0.0,  T: 0.0},},
//                 },
//             },
//         Name1: "chr1",
// 		Input2: []PFasta{
//             PFasta{Name: "chr1",
//                 Seq: []pDna.Float32Base{
//                     {A: 0.2,  C: 0.3,  G: 0.4,  T: 0.1},
//                     {A: 0.25, C: 0.25, G: 0.25, T: 0.25},
//                     {A: 0.2,  C: 0.3,  G: 0.3,  T: 0.2},
//                     {A: 0.1,  C: 0.2,  G: 0.3,  T: 0.4},
//                     {A: 0.6,  C: 0.2,  G: 0.1,  T: 0.1},},
//                 },
//             PFasta{Name: "chr2",
//                 Seq: []pDna.Float32Base{
//                     {A: 0.5,  C: 0.5,  G: 0.0,  T: 0.0},
//                     {A: 0.5,  C: 0.5,  G: 0.0,  T: 0.0},
//                     {A: 0.5,  C: 0.5,  G: 0.0,  T: 0.0},
//                     {A: 0.5,  C: 0.5,  G: 0.0,  T: 0.0},
//                     {A: 0.5,  C: 0.5,  G: 0.0,  T: 0.0},},
//                     },
//                 },
//         Name2: "chr1",
// 		OutName:      "chr1.dist",
// 		DefaultValue: 0.0,
// 		Expected: map[string]wig.Wig{
// 			"chr1": wig.Wig{
//                 StepType:     "fixedStep",
//                 Chrom:        "chr1.dist",
//                 Start:        1,
//                 Step:         1,
//                 Span:         1,
//                 DefaultValue: 0.0,
//                 Values: []float64{0, 0, 0, 0, 0},},},
//         Precision: 0.001,
//     },
// }


// func TestDistTrackMulti(t *testing.T) {
//     var res map[string]wig.Wig
// 	for _, v := range DistTrackMultiTests {
// 		res = map[string]wig.Wig{"chr1": DistTrackMulti(v.Input1, v.Name1, v.Input2, v.Name2, v.OutName, v.DefaultValue)}
// 		if !wig.AllEqual(res, v.Expected, v.Precision) {			t.Errorf("Error: in pFasta. DistTrack valid input test was not as expected.")
// 		}
// 	}
// }
