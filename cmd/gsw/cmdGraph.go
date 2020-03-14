package main

import (
	"bytes"
	"flag"
	"fmt"
	"github.com/vertgenlab/gonomics/axt"
	"github.com/vertgenlab/gonomics/fasta"
	"github.com/vertgenlab/gonomics/sam"
	"github.com/vertgenlab/gonomics/simpleGraph"
	"github.com/vertgenlab/gonomics/vcf"
	"log"
	"os"
	"os/exec"
	"strings"
)

func usage() {
	fmt.Print(
		"\nGSW - genome graph toolkit\n" +
			"\nusage:\n" +
			"\tgsw [options] ref.fa/.gg R1.fastq.gz R2.fastq.gz\n" +
			"options:\n" +
			"\t--align\n" +
			"\t\tGraph-Smith-Waterman: align fastqs to genome graph\n" +
			"\t\t./gsw --align --out align.sam ref.gg R1.fastq.gz R2.fastq.gz\n" +
			"\t--ggTools\n" +
			"\t\tcreate genome graph reference w/ vcf file\n" +
			"\t\t./gsw --vcf input.vcf --out ref.gg ref.fa\n" +
			"\t--view\n" +
			"\t\tvisualize alignment /dev/stdout\n" +
			"\t\t./gsw --view alignToGraph.sam ref.gg\n" +
			"\t--options\n" +
			"\t\tneed more help? See advanced user options:\n" +
			"\t\t./gsw --options [align/ggTools/view]\n\n")
}
func needHelp(cmdName string) {
	var answer string = "" //usage:\n\t//gsw [options] ref.fa/.gg R1.fastq.gz R2.fastq.gz\noptions:\n"
	if strings.Compare(cmdName, "align") == 0 {
		answer += "GSW - genome graph toolkit:\n\n" +
			"Graph-Smith-Waterman:\n\n" + "usage:" +
			"\t./gsw --align [options] ref.gg R1.fastq.gz R2.fastq.gz\n\n" +
			"\t--seed\tdefault: 32\n" +
			"\t\tseedLen kMer for creating hash look up reference genome\n" +
			"\t\t> kMer: slower alignment, better mapping\n" +
			"\t\t< kMer: faster alignment, but higher change of unmapped reads\n" +
			"\t\tbetween 2 and 32\n\n" +
			"\t--step\tdefault: k-1\n" +
			"\t\toffset position of sliding window of hash\n\n" +
			"\t--cpu\tdefault: 4\n" +
			"\t\tnumber of CPUs to use\n\n"
	} else if strings.Compare(cmdName, "ggTools") == 0 {
		answer += "\nggTools: utilities to create, manipulate and operate on genome graphs\n" +
			"\nTo create genome graph reference w/ vcf file:\n\n" +
			"./gsw --ggTools --vcf SNPsIndels.vcf genome.fa\n\n" +
			"usage:\t./gsw --ggTools [options] ref.[.gg/.fa]\n\n" +
			"\t--split\t--out genome_[chr1, chr2, chr3 ...].gg\n" +
			"\t\tgraph reference split by chromosome\n\n" +
			"\t--merge\t./gsw --ggTools merged.sam --merge chr1.sam chr2.sam chr3.sam...\n" +
			"\t\tmerge split by chromosome sam files into one\n" +
			"\t\tfinds best alignment for each read\n\n" +
			"\t--axt\t./gsw --ggTools --axt genomes.axt --out SNPsIndels.vcf ref.fa\n" +
			"\t\tuse axt alignment to create VCF: small SNPs and indels\n\n" +
			"\t--slurm\tbeta: submit GSW command as a slurm job\n\n"
	} else if strings.Compare(cmdName, "view") == 0 {
		answer += "GSW - genome graph toolkit:\n\n" + "visualize alignment\n" + "usage:" +
			"\t./gsw --view [options] --vcf input.vcf ref.[.gg/.fa]\n\n" +
			"\t--out\tdefault: /dev/stdout\n" +
			"\t\tnotes.txt\n"
	} else if strings.Compare(cmdName, "axt") == 0 {
		answer +=
			"\t./gsw  --out SNPsIndels.vcf ref.fa\n\n"
	} else {
		log.Fatalf("Error: Apologies, your command line prompt was not recognized...")
	}
	fmt.Print(answer)
}

//WORK IN PROGRESS
func main() {
	var tagAxt *string = flag.String("axt", "", "axt alignment file")
	var vcfTag *string = flag.String("vcf", "", "vcf file")
	var outTag *string = flag.String("out", "/dev/stdout", "final output, .vcf/.gg/.sam")
	var alignFlag *bool = flag.Bool("align", false, "in.fastq out.sam")
	var threads *int = flag.Int("cpu", 4, "Number of threads or CPUs to use")
	var kMerHash *int = flag.Int("seed", 32, "Seed length used for indexing the reference genome")
	var stepSize *int = flag.Int("step", *kMerHash-1, "step size for building hash")
	var view *string = flag.String("view", "", "visualize sam alignment")
	var moreHelp *string = flag.String("options", "", "advanced user options")
	var splitChr *bool = flag.Bool("split", false, "splits graph output by chromosomes")
	var chrPrefix *string = flag.String("name", "genomeGraph", "basename for .gg file, split by chromosome")

	var ggTools *bool = flag.Bool("ggTools", false, "genome graph tools")
	var mergeSam *string = flag.String("merge", "", "merge split sam files back into one")
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
		}
	}
	//log.Printf("Num of args=%d\n", len(flag.Args()))
	if *slurmScript == true && *ggTools == true {
		slurm()
	} else if *kent == true && *ggTools == true {
		kentUtils(flag.Args())
	} else {
		if *alignFlag == true {
			log.Printf("Reading reference genome...\n")
			ref, chrSize := simpleGraph.Read(flag.Arg(0))
			//user provides single end reads
			if len(flag.Args()) == 2 {
				
				simpleGraph.GswSingleReadWrap(ref, flag.Arg(1), *outTag, *threads, *kMerHash, *stepSize, chrSize)
			} else if len(flag.Args()) == 3 {
				//user provides paired end reads
				simpleGraph.GSWsBatchPair(ref, flag.Arg(1), flag.Arg(2), *outTag, *threads, *kMerHash, chrSize)
			}
		}
		if *ggTools == true {
			if strings.HasSuffix(*tagAxt, ".axt") {
				if *outTag != "" {
					axtFile := axt.Read(*tagAxt)
					fa := fasta.Read(flag.Arg(0))
					axt.AxtVcfToFile(*outTag, axtFile, fa)
				}
			} else if strings.HasSuffix(*vcfTag, ".vcf") {
				vcfs := vcf.Read(*vcfTag)
				if strings.HasSuffix(flag.Arg(0), ".fa") {
					fa := fasta.Read(flag.Arg(0))
					if *splitChr == true {
						ggChr := simpleGraph.SplitGraphChr(fa, vcfs)
						simpleGraph.WriteToGraphSplit(*chrPrefix, ggChr)
					} else {
						gg := simpleGraph.VariantGraph(fa, vcfs)
						simpleGraph.Write(*outTag, gg)
					}
				} else {
					log.Fatalf("Error: Apologies, your command line prompt was not recognized...\nPlease provide both a fasta reference and a VCF file...\n")
				}
			} else if strings.Compare(*mergeSam, "") != 0 {
				sam.ReadFiles(flag.Args(), *mergeSam)
			}
		}
		if strings.HasSuffix(*view, ".sam") {
			var yes, no, numReads int = 0, 0, 0
			log.SetFlags(log.Ldate | log.Ltime)
			if strings.HasSuffix(flag.Arg(0), ".gg") {
				gg, _ := simpleGraph.Read(flag.Arg(0))
				samfile, _ := sam.Read(*view)
				for _, samline := range samfile.Aln {
					log.Printf("%s\n", simpleGraph.ViewGraphAlignment(samline, gg))
					numReads++
					if simpleGraph.CheckAlignment(samline) {
						yes++
					} else {
						no++
					}
				}
				log.Printf("Total number of reads aligned: %d...", numReads)
				log.Printf("Number of reads correctly aligned: %d...\n", yes)
				log.Printf("Number of reads mismapped: %d...\n", no)
			}
		}
	}
}

func slurm() {
	//set basic commands for now:
	slurmJob := "sbatch"
	args := []string{"--mem=32G", "--nodes=1", "--ntasks=1", "--cpus-per-task=8", "--mail-type=END,FAIL", "--mail-user=eric.au@duke.edu"}
	commands := os.Args
	var tmp []string
	for i := 0; i < len(commands); i++ {
		if !strings.Contains(commands[i], "slurm") {
			tmp = append(tmp, commands[i])
		}
	}
	var wrapPrompt string = "--wrap=\""
	wrapPrompt += strings.Join(tmp, " ") + "\""
	args = append(args, wrapPrompt)
	fmt.Printf("The wrap prompt is: %s\n", wrapPrompt)
	cmd := exec.Command(slurmJob, args...)

	//cmd := exec.Command("sbatch", mem, inputMem, nodes, inputNodes, tasks, inputTasks, cpus, numCpus, mailType, email, userEmail, wrap, gswCmd)
	cmdOutput := &bytes.Buffer{}
	cmd.Stdout = cmdOutput
	err := cmd.Run()
	if err != nil {
		os.Stderr.WriteString(err.Error())
	}
	fmt.Print(string(cmdOutput.Bytes()))
}

func kentUtils(command []string) {
	dir := "/Users/bulbasaur/kentUtils/"
	dir += command[0]
	var args []string = command[1:]
	cmd := exec.Command(dir, args...)
	cmdOutput := &bytes.Buffer{}
	cmd.Stdout = cmdOutput
	err := cmd.Run()
	if err != nil {
		os.Stderr.WriteString(err.Error())
	}
	fmt.Print(string(cmdOutput.Bytes()))
	

}
