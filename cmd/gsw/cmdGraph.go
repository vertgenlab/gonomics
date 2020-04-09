package main

import (
	"bytes"
	"flag"
	"fmt"
	"github.com/vertgenlab/gonomics/align"
	"github.com/vertgenlab/gonomics/axt"
	"github.com/vertgenlab/gonomics/chromInfo"
	"github.com/vertgenlab/gonomics/fasta"
	"github.com/vertgenlab/gonomics/sam"
	"github.com/vertgenlab/gonomics/simpleGraph"
	"github.com/vertgenlab/gonomics/vcf"
	"log"
	"os"
	"os/exec"
	"path"
	"path/filepath"
	"strings"
)

func usage() {
	fmt.Print(
		"\nGSW - genome graph toolkit" +
			"\nusage:\n" +
			"\t./gsw [options] ref.fa/.gg\n" +
			"options:\n" +
			"\t--align\n" +
			"\t\tGraph-Smith-Waterman: align fastqs to genome graph\n" +
			"\t\t./gsw --align --out align.sam ref.gg R1.fastq.gz R2.fastq.gz\n" +
			"\t--ggTools\n" +
			"\t\tCreate genome graph reference w/ vcf file\n" +
			"\t\t./gsw --ggTools --vcf input.vcf --out ref.gg ref.fa\n" +
			"\t--view\n" +
			"\t\tVisualize alignment /dev/stdout\n" +
			"\t\t./gsw --view alignToGraph.sam ref.gg\n" +
			"\t--options\n" +
			"\t\tNeed more help? See advanced user options:\n" +
			"\t\t./gsw --options [align/ggTools/view]\n\n")
}
func needHelp(cmdName string) {
	var answer string = "" //usage:\n\t//gsw [options] ref.fa/.gg R1.fastq.gz R2.fastq.gz\noptions:\n"
	if strings.Compare(cmdName, "align") == 0 {
		answer += "GSW - genome graph toolkit:\n\n" +
			"Graph-Smith-Waterman:\n\n" + "usage:" +
			"\t./gsw --align [options] ref.gg R1.fastq.gz R2.fastq.gz\n\n" +
			"\t--seed\tdefault: 32\n" +
			"\t\t\tseedLen kMer for creating hash look up reference genome\n" +
			"\t\t\t> kMer: slower alignment, better mapping\n" +
			"\t\t\t< kMer: faster alignment, but higher change of unmapped reads\n" +
			"\t\t\tbetween 2 and 32\n\n" +
			"\t--step\t\tdefault: 32\n" +
			"\t\t\toffset position of sliding window of hash\n\n" +
			"\t--score\t\thumanChimp:\n\t\t" + printMatrix(align.HumanChimpTwoScoreMatrix) + "\n" +
			"\t\t\thoxD55:\n\t\t" + printMatrix(align.HoxD55ScoreMatrix) + "\n" +
			"\t\t\tmouseRat:\n\t\t" + printMatrix(align.MouseRatScoreMatrix) + "\n" +
			"\t--cpus\t\tdefault: 4\n" +
			"\t\t\tnumber of CPUs to use\n\n" +
			"\t--liftover\tprovide a chrom sizes file\n" +
			"\t\t\tconvert alignment coordinates to linear reference in sam format\n\n"

	} else if strings.Compare(cmdName, "ggTools") == 0 {
		answer += "\nggTools: utilities to create, manipulate and operate on genome graphs\n" +
			"\nTo create genome graph reference w/ vcf file:\n" +
			"./gsw --ggTools --vcf SNPsIndels.vcf genome.fa\n\n" +
			"usage:\t./gsw --ggTools [options] ref.[.gg/.fa]\n\n" +
			"\t--vcf\tSNPsIndels.vcf --out ref.gg ref.fa\n\n" +
			"\t--split\t--out genome_[chr1, chr2, chr3 ...].gg\n" +
			"\t\tgraph reference split by chromosome\n\n" +
			"\t--merge\t--out merge.sam chr1.sam chr2.sam chr3.sam ...\n" +
			"\t\tmerge split by chromosome sam files into one\n" +
			"\t\tfinds best alignment for each read\n\n" +
			"\t--axt\tgenomes.axt --out SNPsIndels.vcf ref.fa\n" +
			"\t\tuse axt alignment to create VCF: small SNPs and indels\n\n" +
			"\t--fa\t--out graph.fa ref.gg\n" +
			"\t\tgenome graph to fasta format\n\n" +
			"\t--slurm\tbeta: submit GSW command as a slurm job\n" +
			"\t\tdefault settings are: --mem=32G, --ntasks=1, --cpus-per-task=8\n\n"
		//"\t--kent\tkentUtils - ucsc\n\n"
	} else if strings.Compare(cmdName, "view") == 0 {
		answer += "GSW - genome graph toolkit:\n\n" + "Visualize Alignment\n\n" + "usage:" +
			"\t./gsw --view [options] graph.sam ref.[.gg/.fa]\n\n" +
			"\t--out\tdefault: /dev/stdout\n" +
			"\t\tnotes.txt\n\n" //+
		//"\tscores\tprint score matrix options\n\n"
	} else if strings.Compare(cmdName, "axt") == 0 {
		answer +=
			"\t./gsw  --out SNPsIndels.vcf ref.fa\n\n"
	} else {
		errorMessage()
		//log.Fatalf("Error: Apologies, your command prompt was not recognized...\n\t\t\t\t\t\t\t\t-xoxo GG")
	}
	fmt.Print(answer)
}

//WORK IN PROGRESS
func main() {
	var tagAxt *string = flag.String("axt", "", "axt alignment file")
	var vcfTag *string = flag.String("vcf", "", "vcf file")
	var outTag *string = flag.String("out", "/dev/stdout", "final output, .vcf/.gg/.sam")
	var alignFlag *bool = flag.Bool("align", false, "in.fastq out.sam")
	var threads *int = flag.Int("cpus", 4, "Number of threads or CPUs to use")
	var kMerHash *int = flag.Int("seed", 32, "Seed length used for indexing the reference genome")
	var stepSize *int = flag.Int("step", 32, "step size for building hash")
	var view *string = flag.String("view", "", "visualize sam alignment")
	var moreHelp *string = flag.String("options", "", "advanced user options")
	var splitChr *bool = flag.Bool("split", false, "splits graph output by chromosomes")
	var chrPrefix *string = flag.String("name", "genomeGraph", "basename for .gg file, split by chromosome")
	var ggTools *bool = flag.Bool("ggTools", false, "genome graph tools")
	var mergeSam *bool = flag.Bool("merge", false, "merge split sam files back into one")
	var score *string = flag.String("score", "humanChimp", "choose a scoring matrix")
	var fa *bool = flag.Bool("fa", false, "converts graph back to fa")
	var liftover *string = flag.String("liftover", "", "liftover to linear reference sam file")

	var slurmScript *bool = flag.Bool("slurm", false, "submit gsw command as a slurm job")
	var kent *bool = flag.Bool("kent", false, "run a kentUtils through GSW")

	flag.Usage = usage
	log.SetFlags(log.Ldate | log.Ltime)
	flag.Parse()
	if *moreHelp != "" {
		needHelp(*moreHelp)
	} else {
		if len(flag.Args()) == 0 {
			flag.Usage()
			errorMessage()
		}
	}
	if *slurmScript && *ggTools && !*splitChr {
		slurm()
	} else if *ggTools && *kent {
		kentUtils(flag.Args())
	} else {
		if *alignFlag {
			log.Printf("Reading reference genome...\n")
			ref := simpleGraph.Read(flag.Arg(0))
			//header.Text = append(header.Text, fmt.Sprintf("@PG\tID:GSW\tPN:ggTools\tVN:1130\tCL:%s", strings.Join(os.Args, " ")))
			//user provides single end reads
			if len(flag.Args()) == 2 {
				simpleGraph.GswSingleReadWrap(ref, flag.Arg(1), *outTag, *threads, *kMerHash, *stepSize, selectScoreMatrix(*score), nil)
			} else if len(flag.Args()) == 3 {
				//user provides paired end reads
				if strings.HasSuffix(*liftover, ".sizes") {
					chrSize := chromInfo.ReadToSlice(*liftover)
					header := sam.ChromInfoSamHeader(chrSize)
					header.Text = append(header.Text, fmt.Sprintf("@PG\tID:GSW\tPN:ggTools\tVN:1130\tCL:%s", strings.Join(os.Args, " ")))
					simpleGraph.WrapGirafLiftoverToSam(ref, flag.Arg(1), flag.Arg(2), *outTag, *threads, *kMerHash, *stepSize, selectScoreMatrix(*score), header)
				} else {
					simpleGraph.GswToGirafPair(ref, flag.Arg(1), flag.Arg(2), *outTag, *threads, *kMerHash, *stepSize, selectScoreMatrix(*score))
				}
			}
		}
		if *ggTools && strings.HasSuffix(*tagAxt, ".axt") {
			if *outTag != "" {
				axtFile := axt.Read(*tagAxt)
				fa := fasta.Read(flag.Arg(0))
				axt.AxtVcfToFile(*outTag, axtFile, fa)
			}
		}
		if *ggTools && (strings.HasSuffix(*vcfTag, ".vcf") || strings.HasSuffix(*vcfTag, ".gz")) {
			vcfs := vcf.Read(*vcfTag)
			if strings.HasSuffix(flag.Arg(0), ".fa") {
				fa := fasta.Read(flag.Arg(0))
				if *splitChr {
					log.Printf("VCF to graph, split into chromosomes...\n")
					ggChr := simpleGraph.SplitGraphChr(fa, vcfs)
					if *slurmScript && strings.Contains(flag.Arg(1), ".fastq") {
						var readOne string = filepath.Base(strings.TrimSuffix(flag.Arg(1), path.Ext(flag.Arg(1))))
						gswCommand := fmt.Sprintf(" --wrap=\"./gsw --align --k 16 --t 8 --out %s_", readOne)
						var currCmd string
						var echo string
						for chr := range ggChr {
							log.Printf("Writing graph %s to file...\n", chr)
							splitRefName := fmt.Sprintf("%s_%s.gg", chr, *outTag)
							simpleGraph.Write(splitRefName, ggChr[chr])

							currCmd = gswCommand + fmt.Sprintf("to_%s_%s.sam %s %s", chr, *outTag, splitRefName, flag.Arg(1))
							if len(flag.Args()) == 3 {
								currCmd += fmt.Sprintf(" %s", flag.Arg(2))
							}
							currCmd += "\""
							args := []string{"--mem=16G", "--nodes=1", "--ntasks=1", "--cpus-per-task=8", "--mail-type=END,FAIL", "--mail-user=eric.au@duke.edu"}
							args = append(args, currCmd)
							log.Printf("\n\nSlurm job submission:\n")
							echo = "sbatch " + strings.Join(args, " ")
							log.Printf("\n\n%s\n\n", echo)
							cmd := exec.Command("sbatch", args...)
							cmdOutput := &bytes.Buffer{}
							cmd.Stdout = cmdOutput
							err := cmd.Run()
							if err != nil {
								os.Stderr.WriteString(err.Error())
							}
							fmt.Print(string(cmdOutput.Bytes()))
						}
					} else {
						simpleGraph.WriteToGraphSplit(*chrPrefix, ggChr)
					}
				} else {
					gg := simpleGraph.VariantGraph(fa, vcfs)
					simpleGraph.Write(*outTag, gg)
				}
			} else {
				errorMessage()
			}
		}
		if *ggTools && *fa {
			gg := simpleGraph.Read(flag.Arg(0))

			f := simpleGraph.GraphToFa(gg)
			fasta.Write(*outTag, f)

		}
		if *ggTools && *mergeSam {
			sam.ReadFiles(flag.Args(), *outTag)
		}
		if strings.HasSuffix(*view, ".sam") {
			log.SetFlags(log.Ldate | log.Ltime)
			if strings.HasSuffix(flag.Arg(0), ".gg") || strings.HasSuffix(flag.Arg(0), ".fa") {
				gg := simpleGraph.Read(flag.Arg(0))
				samfile, _ := sam.Read(*view)
				for _, samline := range samfile.Aln {
					log.Printf("%s\n", simpleGraph.ViewGraphAlignment(samline, gg))
				}

			}
		}
	}
}

func errorMessage() {
	log.Fatalf("Error: Apologies, your command prompt was not recognized...\n\n-xoxo GG\n")
}

func selectScoreMatrix(score string) [][]int64 {
	var scoreMatrix [][]int64
	switch {
	case strings.Contains(score, "humanChimp"):
		scoreMatrix = align.HumanChimpTwoScoreMatrix
	case strings.Contains(score, "hoxD55"):
		scoreMatrix = align.HoxD55ScoreMatrix
	case strings.Contains(score, "mouseRat"):
		scoreMatrix = align.MouseRatScoreMatrix
	case strings.Contains(score, "general"):
		scoreMatrix = align.DefaultScoreMatrix
	default:
	}
	return scoreMatrix
}

func printMatrix(m [][]int64) string {
	var message string = ""
	message += fmt.Sprintf("\t %d\t%d\t%d\t%d\n\t\t\t%d\t%d\t%d\t%d\n\t\t\t%d\t%d\t%d\t%d\n\t\t\t%d\t%d\t%d\t %d\n", m[0][0], m[0][1], m[0][2], m[0][3], m[1][0], m[1][1], m[1][2], m[1][3], m[2][0], m[2][1], m[2][2], m[2][3], m[3][0], m[3][1], m[3][2], m[3][3])
	return message
}

//TODO: Will remove to a personal script
func slurm() {
	//set basic commands for now:
	slurmJob := "sbatch"
	args := []string{"--mem=32G", "--nodes=1", "--ntasks=1", "--cpus-per-task=8"}
	//, "--mail-type=END,FAIL", "--mail-user=eric.au@duke.edu"
	var gswCommand []string
	var wrapPrompt string = "--wrap=\""

	for _, cmds := range os.Args {
		if !strings.Contains(cmds, "slurm") {
			gswCommand = append(gswCommand, cmds)
		}
	}
	wrapPrompt += strings.Join(gswCommand, " ") + "\""
	args = append(args, wrapPrompt)
	echo := slurmJob + " " + strings.Join(args, " ")
	log.Printf("\n\nSlurm job submission:\n\n%s\n\n", echo)
	cmd := exec.Command(slurmJob, args...)
	cmdOutput := &bytes.Buffer{}
	cmd.Stdout = cmdOutput
	err := cmd.Run()
	if err != nil {
		os.Stderr.WriteString(err.Error() + "\n")
	}
	fmt.Print(string(cmdOutput.Bytes()) + "\n")
}

//TODO: Will remove to a personal script
func kentUtils(command []string) {
	dir := "/Users/edotau/kentUtils/"
	if len(command) == 0 {
		cmd := exec.Command("ls", dir)
		cmdOutput := &bytes.Buffer{}
		cmd.Stdout = cmdOutput
		err := cmd.Run()
		if err != nil {
			os.Stderr.WriteString(err.Error())
		}
		fmt.Print(string(cmdOutput.Bytes()))
	} else {
		dir += command[0]
		var args []string = command[1:]
		cmd := exec.Command(dir, args...)
		cmdOutput := &bytes.Buffer{}
		cmd.Stdout = cmdOutput
		err := cmd.Run()
		if err != nil {
			os.Stderr.WriteString(err.Error() + "\n")
		}
		fmt.Print(string(cmdOutput.Bytes()) + "\n")
	}
}
