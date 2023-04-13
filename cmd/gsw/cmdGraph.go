package main

import (
	"flag"
	"fmt"
	"os"
)

// Genome Graph Gonomics
func usage() {
	fmt.Print(
		"GSW - Graph Smith Waterman:\n\nGenome Graph Software\n" +
			"Author: Eric Au\n\tCraig Lowe\n\n" +
			"Vertebrate Genetics Laboratory: http://www.vertgenlab.org\n" +
			"Source code: https://github.com/vertgenlab/gonomics\n" +
			"Documents: https://github.com/vertgenlab/vglDocumentation\n\n" +
			"Version: 0.1.0\n\n" +
			"Usage:\n" +
			"  gsw [options]\n\n" +
			"Options:\n")
	alignUsage()
	ggToolsUsage()
	viewUsage()
	helpMessage()
	flagsPrint()
}

func main() {
	flag.Usage = usage
	if len(os.Args) < 2 {
		usage()
	} else {
		switch os.Args[1] {
		case "align":
			if len(os.Args) == 2 { //2 ./gsw align, this catch makes it equivalent to gsw align -h includes the name first command
				alignUsage()
				extendedAlignUsage()
			} else {
				RunAlignExe()
			}
		case "ggtools":
			RunGgTools()
		case "view":
			RunViewExe()
		case "help":
			extendHelpMsg.Parse(os.Args[2:])
			moreHelp(os.Args[2])
		default:
			errorMessage()
		}
	}
}

/*
func main() {
	if err := root(os.Args[1:]); err != nil {
		fmt.Println(err)
		os.Exit(1)
	}
}*/

/*
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
}*/
