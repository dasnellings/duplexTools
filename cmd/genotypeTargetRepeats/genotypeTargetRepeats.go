package main

import (
	"errors"
	"flag"
	"fmt"
	"github.com/dasnellings/MCS_MS/fai"
	"github.com/dasnellings/MCS_MS/gmm"
	"github.com/dasnellings/MCS_MS/realign"
	"github.com/guptarohit/asciigraph"
	"github.com/vertgenlab/gonomics/bed"
	"github.com/vertgenlab/gonomics/cigar"
	"github.com/vertgenlab/gonomics/dna"
	"github.com/vertgenlab/gonomics/exception"
	"github.com/vertgenlab/gonomics/fasta"
	"github.com/vertgenlab/gonomics/fileio"
	"github.com/vertgenlab/gonomics/sam"
	"github.com/vertgenlab/gonomics/vcf"
	"golang.org/x/exp/slices"
	"io"
	"log"
	"math"
	"os"
	"path"
	"runtime"
	"runtime/pprof"
	"sort"
	"strconv"
	"strings"
)

var debug int = 0

func usage() {
	fmt.Print(
		"genotypeTargetRepeats - Output a VCF of genotypes of targeted short simple repeats.\n\n" +
			"options:\n")
	flag.PrintDefaults()
}

// inputFiles is a custom type that gets filled by flag.Parse()
type inputFiles []string

// String to satisfy flag.Value interface
func (i *inputFiles) String() string {
	return strings.Join(*i, " ")
}

// Set to satisfy flag.Value interface
func (i *inputFiles) Set(value string) error {
	*i = append(*i, value)
	return nil
}

func main() {
	var inputs inputFiles
	flag.Var(&inputs, "i", "Input BAM file with alignments. Must be sorted and indexed. Can be declared more than once")
	var ref *string = flag.String("r", "", "Reference genome. Must be the same reference used for generating the BAM file.")
	var targets *string = flag.String("t", "", "BED file of targeted repeats. The 4th column must be the sequence of one repeat unit (e.g. CA for a CACACACA repeat), or 'RepeatLen'x'RepeatSeq' (e.g. 10xCA).")
	var output *string = flag.String("o", "stdout", "Output VCF file.")
	var lenOut *string = flag.String("lenOut", "", "Output a bed file with additional columns for determined read lengths for each sample.")
	var bamOut *string = flag.String("bamOutPfx", "", "Output a BAM file with realigned reads. Only outputs reads that inform called genotypes. File will be named 'bamOutPfx'_'originalFilename'.")
	var targetPadding *int = flag.Int("tPad", 50, "Add INT bases of padding to either end of regions in targets file for selecting reads for realignment.")
	var minFlankOverlap *int = flag.Int("minFlank", 4, "A minimum of INT bases must be mapped on either side of the repeat to be considered an enclosing read.")
	var minMapQ *int = flag.Int("minMapQ", -1, "Minimum mapping quality (before realignment) to be considered for genotyping. Set to -1 for no filter.")
	var allowDups *bool = flag.Bool("allowDups", false, "Do not remove duplicate reads when genotyping.")
	var debugVal *int = flag.Int("debug", 0, "Set to 1 or greater for debug prints.")
	var minReads *int = flag.Int("minReads", 5, "Minimum total enclosing reads for genotyping.")
	var alignerThreads *int = flag.Int("alnThreads", 1, "Number of alignment threads.")
	cpuprofile := flag.String("cpuprofile", "", "write cpu profile to `file`")
	memprofile := flag.String("memprofile", "", "write memory profile to `file`")
	flag.Parse()
	flag.Usage = usage

	if *cpuprofile != "" {
		f, err := os.Create(*cpuprofile)
		if err != nil {
			log.Fatal("could not create CPU profile: ", err)
		}
		defer f.Close()
		if err := pprof.StartCPUProfile(f); err != nil {
			log.Fatal("could not start CPU profile: ", err)
		}
		defer pprof.StopCPUProfile()
	}

	if len(inputs) == 0 || *ref == "" {
		usage()
		log.Fatalln("ERROR: must input a VCF file with -i")
	}

	debug = *debugVal

	if *minMapQ > math.MaxUint8 {
		log.Fatalf("minMapQ out of range. max: %d\n", math.MaxUint8)
	}

	genotypeTargetRepeats(inputs, *ref, *targets, *output, *bamOut, *lenOut, *targetPadding, *minFlankOverlap, *minMapQ, *minReads, !*allowDups, *alignerThreads)

	if *memprofile != "" {
		f, err := os.Create(*memprofile)
		if err != nil {
			log.Fatal("could not create memory profile: ", err)
		}
		defer f.Close() // error handling omitted for example
		runtime.GC()    // get up-to-date statistics
		if err := pprof.WriteHeapProfile(f); err != nil {
			log.Fatal("could not write memory profile: ", err)
		}
	}
}

func genotypeTargetRepeats(inputFiles []string, refFile, targetsFile, outputFile, bamOutPfx, lenOutFile string, targetPadding, minFlankOverlap, minMapQ, minReads int, removeDups bool, alignerThreads int) {
	var err error
	var ref *fasta.Seeker
	var lenOut *fileio.EasyWriter
	targets := bed.Read(targetsFile)
	vcfOut := fileio.EasyCreate(outputFile)
	defer cleanup(vcfOut)
	vcfHeader := generateVcfHeader(strings.Join(inputFiles, "\t"), refFile)
	vcf.NewWriteHeader(vcfOut, vcfHeader)

	// get bam reader for each file
	br := make([]*sam.BamReader, len(inputFiles))
	headers := make([]sam.Header, len(inputFiles))
	bamIdxs := make([]sam.Bai, len(inputFiles))
	for i := range inputFiles {
		br[i], headers[i] = sam.OpenBam(inputFiles[i])
		defer cleanup(br[i])
		if _, err = os.Stat(inputFiles[i] + ".bai"); !errors.Is(err, os.ErrNotExist) {
			bamIdxs[i] = sam.ReadBai(inputFiles[i] + ".bai")
		} else {
			bamIdxs[i] = sam.ReadBai(strings.TrimSuffix(inputFiles[i], ".bam") + ".bai")
		}
	}

	bamOutHandle := make([]io.WriteCloser, len(inputFiles))
	bamOut := make([]*sam.BamWriter, len(inputFiles))
	if bamOutPfx != "" {
		for i := range inputFiles {
			words := strings.Split(inputFiles[i], "/")
			words[len(words)-1] = bamOutPfx + "_" + words[len(words)-1]
			bamOutHandle[i] = fileio.EasyCreate(words[len(words)-1])
			bamOut[i] = sam.NewBamWriter(bamOutHandle[i], headers[i])
			defer cleanup(bamOutHandle[i])
			defer cleanup(bamOut[i])
		}
	}

	if lenOutFile != "" {
		lenOut = fileio.EasyCreate(lenOutFile)
		fmt.Fprintf(lenOut, "#CHROM\tSTART\tEND\tREPEAT\t%s\n", strings.Join(inputFiles, "\t"))
		defer cleanup(lenOut)
	}

	enclosingReads := make([][]*sam.Sam, len(inputFiles)) // first index is sample
	observedLengths := make([][]int, len(inputFiles))     // first index is sample
	var currVcf vcf.Vcf
	alignerInput := make(chan sam.Sam, 100)
	alignerOutput := make(chan sam.Sam, 100)
	for j := 0; j < alignerThreads; j++ {
		ref = fasta.NewSeeker(refFile, "")
		defer cleanup(ref)
		go realign.RealignIndels(alignerInput, alignerOutput, ref)
	}

	mm := make([]*gmm.MixtureModel, len(inputFiles))
	for i := 0; i < len(inputFiles); i++ {
		mm[i] = new(gmm.MixtureModel)
	}
	tmpMm := new(gmm.MixtureModel)

	gaussians := make([][]float64, 2)
	var floatSlice []float64
	var converged, anyConverged, passingVariant bool
	for _, region := range targets {
		anyConverged = false
		for i := range inputFiles {
			enclosingReads[i], observedLengths[i] = getLenghtDist(enclosingReads[i], targetPadding, minMapQ, minFlankOverlap, removeDups, bamIdxs[i], region, br[i], bamOut[i], alignerInput, alignerOutput)
			if bamOutPfx != "" {
				for j := range enclosingReads[i] {
					sam.WriteToBamFileHandle(bamOut[i], *enclosingReads[i][j], 0)
				}
			}
			slices.Sort(observedLengths[i])

			converged, tmpMm, mm[i] = runMixtureModel(observedLengths[i], tmpMm, mm[i], &floatSlice)
			if converged {
				anyConverged = true
			}
		}

		if !anyConverged {
			continue
		}

		if lenOut != nil {
			fmt.Fprintf(lenOut, "%s%s\n", bed.ToString(region, 4), printLengths(observedLengths))
		}

		if debug > 0 {
			plot(observedLengths, minReads, mm, gaussians)
		}

		currVcf, passingVariant = callGenotypes(ref, region, minReads, enclosingReads, observedLengths, mm)
		if passingVariant {
			vcf.WriteVcf(vcfOut, currVcf)
		}
	}
	close(alignerInput)
	close(alignerOutput)
}

func callGenotypes(ref *fasta.Seeker, region bed.Bed, minReads int, enclosingReads [][]*sam.Sam, observedLengths [][]int, mm []*gmm.MixtureModel) (vcf.Vcf, bool) {
	var ans vcf.Vcf
	repeatUnitLen, refNumRepeats := parseRepeatSeq(region.Name)
	refRepeatLen := refNumRepeats * len(repeatUnitLen)
	ans.Chr = region.Chrom
	ans.Pos = region.ChromStart
	refSeq, err := fasta.SeekByName(ref, region.Chrom, region.ChromStart, region.ChromEnd)
	exception.PanicOnErr(err)
	dna.AllToUpper(refSeq)
	ans.Ref = dna.BasesToString(refSeq)
	ans.Ref = "*" // TODO Remove
	//if len(ans.Ref) != refRepeatLen {
	//	log.Panicf("ERROR: %s ref seq is \n%s\n the length of %d does not match expected %d from bed file.", region, ans.Ref[1:], len(ans.Ref), refRepeatLen)
	//}

	ans.Id = region.Name

	/*
		altLens := make([]int, 2)
		var refLenDiff int
		for i, l := range mm[0].Means {
			altLens[i] = int(math.Round(l))
			refLenDiff = refRepeatLen - altLens[i]
			for _, alts := range ans.Alt {
				if len(alts) == altLens[i] {
					refLenDiff = 0 // to engage break below
				}
			}
			if refLenDiff == 0 {
				continue
			}
			ans.Alt = append(ans.Alt, ans.Ref[0:len(ans.Ref)-refLenDiff-1])
		}
	*/
	ans.Alt = append(ans.Alt, "*")
	ans.Filter = "."
	ans.Id = region.Name
	ans.Format = []string{"GT", "DP", "MU", "SD", "WT", "LL"}
	ans.Samples = make([]vcf.Sample, len(mm))
	for i := range ans.Samples {
		ans.Samples[i].FormatData = make([]string, 6)
		ans.Samples[i].FormatData[1] = fmt.Sprintf("%d", len(observedLengths[i]))
		if mm[i].LogLikelihood == math.MaxFloat64 {
			ans.Samples[i].FormatData[2] = "."
			ans.Samples[i].FormatData[3] = "."
			ans.Samples[i].FormatData[4] = "."
			ans.Samples[i].FormatData[5] = "."
			continue
		}

		if mm[i].Means[0] < mm[i].Means[1] {
			ans.Samples[i].FormatData[2] = fmt.Sprintf("%.1f,%.1f", mm[i].Means[0], mm[i].Means[1])
			ans.Samples[i].FormatData[3] = fmt.Sprintf("%.1f,%.1f", mm[i].Stdev[0], mm[i].Stdev[1])
			ans.Samples[i].FormatData[4] = fmt.Sprintf("%.1f,%.1f", mm[i].Weights[0], mm[i].Weights[1])
		} else {
			ans.Samples[i].FormatData[2] = fmt.Sprintf("%.1f,%.1f", mm[i].Means[1], mm[i].Means[0])
			ans.Samples[i].FormatData[3] = fmt.Sprintf("%.1f,%.1f", mm[i].Stdev[1], mm[i].Stdev[0])
			ans.Samples[i].FormatData[4] = fmt.Sprintf("%.1f,%.1f", mm[i].Weights[1], mm[i].Weights[0])
		}
		ans.Samples[i].FormatData[5] = fmt.Sprintf("%.1g", mm[i].LogLikelihood)
	}

	info := new(strings.Builder)
	info.WriteString(fmt.Sprintf("RefLength=%d", refRepeatLen))
	/*
		info.WriteString(";Means=")
		for i, j := range mm[0].Means {
			if i > 0 {
				info.WriteByte(',')
			}
			info.WriteString(fmt.Sprintf("%.1f", j))
		}
		info.WriteString(";Stdev=")
		for i, j := range mm[0].Stdev {
			if i > 0 {
				info.WriteByte(',')
			}
			info.WriteString(fmt.Sprintf("%.1f", j))
		}
		info.WriteString(";Weights=")
		for i, j := range mm[0].Weights {
			if i > 0 {
				info.WriteByte(',')
			}
			info.WriteString(fmt.Sprintf("%.1f", j))
		}
	*/
	ans.Info = info.String()

	return ans, true
}

func getLenghtDist(enclosingReads []*sam.Sam, targetPadding, minMapQ, minFlankOverlap int, removeDups bool, bamIdx sam.Bai, region bed.Bed, br *sam.BamReader, bamOut *sam.BamWriter, alignerInput chan<- sam.Sam, alignerOutput <-chan sam.Sam) ([]*sam.Sam, []int) {
	var start, end int
	var reads []sam.Sam
	enclosingReads = resetEnclosingReads(enclosingReads, len(reads)) // starts at len == 0, cap >= len(reads)

	// STEP 1: Find reads with initial alignment close to target as candidates for local realignment
	start = region.ChromStart - targetPadding
	end = region.ChromEnd + targetPadding
	if start < 0 {
		start = 0
	}
	reads = sam.SeekBamRegion(br, bamIdx, region.Chrom, uint32(start), uint32(end))
	if len(reads) == 0 {
		return enclosingReads, nil
	}

	// STEP 2: Realign reads to target region
	realignReads(reads, minMapQ, alignerInput, alignerOutput) // read order in slice may change

	// STEP 3: Determine which realigned reads overlap targets with the minimum flanking overlap
	for i := range reads {
		if minMapQ != -1 && reads[i].MapQ < uint8(minMapQ) {
			continue
		}
		if sam.IsUnmapped(reads[i]) {
			continue
		}
		if reads[i].GetChromStart() <= region.ChromStart-minFlankOverlap && reads[i].GetChromEnd() >= region.ChromEnd+minFlankOverlap {
			enclosingReads = append(enclosingReads, &reads[i])
		}
	}

	// STEP 4: Sort enclosing reads by position
	sort.Slice(enclosingReads, func(i, j int) bool {
		if enclosingReads[i].GetChromStart() < enclosingReads[j].GetChromStart() {
			return true
		}
		if enclosingReads[i].GetChromEnd() < enclosingReads[j].GetChromEnd() {
			return true
		}
		return true
	})

	// STEP 5: Remove duplicates
	if removeDups {
		enclosingReads = dedup(enclosingReads)
	}

	// STEP 6: Genotype repeats
	observedLengths := make([]int, len(enclosingReads))
	repeatSeq, _ := parseRepeatSeq(region.Name)
	for i := range enclosingReads {
		observedLengths[i] = calcRepeatLength(enclosingReads[i], region.ChromStart, region.ChromEnd, repeatSeq)
		if debug > 2 {
			fmt.Fprintln(os.Stderr, enclosingReads[i].QName, observedLengths[i], "start:", enclosingReads[i].Pos)
		}
	}
	return enclosingReads, observedLengths
}

func calcRepeatLength(read *sam.Sam, regionStart, regionEnd int, repeatSeq []dna.Base) int {
	var readIdx, refIdx, i int
	refIdx = int(read.Pos)

	// get to start of region
	for i = range read.Cigar {
		if cigar.ConsumesReference(read.Cigar[i].Op) {
			refIdx += read.Cigar[i].RunLength
		}
		if cigar.ConsumesQuery(read.Cigar[i].Op) {
			readIdx += read.Cigar[i].RunLength
		}
		if refIdx >= regionStart {
			break
		}
	}
	if refIdx > regionStart {
		if cigar.ConsumesQuery(read.Cigar[i].Op) {
			readIdx -= refIdx - regionStart
		}
		refIdx -= refIdx - regionStart
	}
	readIdx++

	var repeatIdx int
	for repeatIdx = range repeatSeq {
		if read.Seq[readIdx] == repeatSeq[repeatIdx] {
			break
		}
	}

	// move backwards to look for misaligned repeat sequence
	for read.Seq[readIdx] == repeatSeq[repeatIdx] {
		repeatIdx--
		readIdx--
		refIdx--
		if repeatIdx == -1 {
			repeatIdx = len(repeatSeq) - 1
		}
		if readIdx == -1 {
			break
		}
	}
	repeatIdx++
	if repeatIdx == len(repeatSeq) {
		repeatIdx = 0
	}
	readIdx++
	refIdx++
	// move forwards to calc repeat length
	var observedLength, maxLength int
	for refIdx < regionEnd && readIdx < len(read.Seq) {
		// move through repeat until mismatch
		for read.Seq[readIdx] == repeatSeq[repeatIdx] {
			observedLength++
			repeatIdx++
			readIdx++
			refIdx++
			if repeatIdx == len(repeatSeq) {
				repeatIdx = 0
			}
			if readIdx == len(read.Seq) {
				break
			}
		}
		if observedLength > maxLength {
			maxLength = observedLength
			observedLength = 0
		}
		// move forward until you get a base matching the repeat
		for readIdx < len(read.Seq) && read.Seq[readIdx] != repeatSeq[repeatIdx] {
			for repeatIdx = 0; repeatIdx < len(repeatSeq); repeatIdx++ {
				if read.Seq[readIdx] == repeatSeq[repeatIdx] {
					break
				}
			}
			if repeatIdx == len(repeatSeq) { // current read base does not match any base in repeat sequence
				repeatIdx = 0
				readIdx++
				refIdx++
			}
		}
	}
	return maxLength // TODO divide by repeat unit length???
}

func parseRepeatSeq(s string) ([]dna.Base, int) {
	var words []string
	if strings.Contains(s, "x") {
		words = strings.Split(s, "x")
	}
	num, err := strconv.Atoi(words[0])
	exception.PanicOnErr(err)
	return dna.StringToBases(words[1]), num
}

func dedup(reads []*sam.Sam) []*sam.Sam {
	for i := 1; i < len(reads); i++ {
		if reads[i].GetChromStart() == reads[i-1].GetChromStart() && reads[i].GetChromEnd() == reads[i-1].GetChromEnd() {
			slices.Delete(reads, i, i+1)
		}
	}
	return reads
}

// read order may change
func realignReads(reads []sam.Sam, minMapQ int, alignerInput chan<- sam.Sam, alignerOutput <-chan sam.Sam) {
	var readsSkipped, readsReceived int

	// count how many reads will be skipped over for realignment
	for i := range reads {
		if minMapQ != -1 && reads[i].MapQ < uint8(minMapQ) {
			readsSkipped++
		}
	}

	// start streaming reads to aligner
	go sendReads(reads, minMapQ, alignerInput)

	// start receiving aligned reads
	for read := range alignerOutput {
		reads[readsReceived] = read
		readsReceived++

		// break when all reads sent for alignment have been received
		if readsReceived+readsSkipped == len(reads) {
			reads = reads[0:readsReceived]
			break
		}
	}
}

func sendReads(reads []sam.Sam, minMapQ int, alignerInput chan<- sam.Sam) {
	for i := range reads {
		if minMapQ != -1 && reads[i].MapQ < uint8(minMapQ) {
			continue
		}
		alignerInput <- reads[i]
	}
}

func resetEnclosingReads(s []*sam.Sam, len int) []*sam.Sam {
	if cap(s) >= len {
		for i := range s {
			s[i] = nil
		}
		s = s[:0]
	} else {
		s = make([]*sam.Sam, 0, len)
	}
	return s
}

func generateVcfHeader(samples string, referenceFile string) vcf.Header {
	var header vcf.Header
	header.Text = append(header.Text, "##fileformat=VCFv4.2")
	header.Text = append(header.Text, fmt.Sprintf("##reference=%s", path.Clean(referenceFile)))
	header.Text = append(header.Text, strings.TrimSuffix(fai.IndexToVcfHeader(fai.ReadIndex(referenceFile+".fai")), "\n"))
	header.Text = append(header.Text, "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">")
	header.Text = append(header.Text, "##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"Total Read Depth\">")
	header.Text = append(header.Text, "##FORMAT=<ID=MU,Number=2,Type=Float,Description=\"Mean repeat length of each allele determined by gaussian mixture modelling.\">")
	header.Text = append(header.Text, "##FORMAT=<ID=SD,Number=2,Type=Float,Description=\"Standard deviation of the repeat length of each allele determined by gaussian mixture modelling.\">")
	header.Text = append(header.Text, "##FORMAT=<ID=WT,Number=2,Type=Float,Description=\"Weight assigned to each allele (rough estimate of allele frequency) determined by gaussian mixture modelling.\">")
	header.Text = append(header.Text, "##FORMAT=<ID=LL,Number=1,Type=Float,Description=\"Negative log likelihood of gaussian mixture model.\">")
	header.Text = append(header.Text, "##INFO=<ID=RefLength,Number=1,Type=Integer,Description=\"Length in bp of the repeat in the reference genome.\">")
	header.Text = append(header.Text, fmt.Sprintf("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t%s", strings.Replace(samples, ".bam", "", -1)))
	return header
}

func plot(observedLengths [][]int, minReads int, mm []*gmm.MixtureModel, gaussians [][]float64) {
	readsPerSample := make([]int, len(observedLengths))
	p := make([][]float64, len(observedLengths))
	for i := range observedLengths {
		p[i] = make([]float64, 100)
		for j := range observedLengths[i] {
			p[i][observedLengths[i][j]]++
			readsPerSample[i]++
		}
	}
	if len(observedLengths) == 1 && readsPerSample[0] < minReads {
		return
	}

	for i := range p {
		if readsPerSample[i] < minReads {
			continue
		}
		//if i != 0 {
		//	continue
		//}
		fmt.Println(asciigraph.Plot(p[i], asciigraph.Height(5), asciigraph.Precision(0), asciigraph.SeriesColors(asciigraph.AnsiColor(i))))

		gaussians[0] = gaussianHist(mm[i].Weights[0], mm[i].Means[0], mm[i].Stdev[0])
		gaussians[1] = gaussianHist(mm[i].Weights[1], mm[i].Means[1], mm[i].Stdev[1])

		fmt.Println(asciigraph.PlotMany(gaussians, asciigraph.Precision(0), asciigraph.SeriesColors(
			asciigraph.Red,
			asciigraph.Yellow,
			asciigraph.Green,
			asciigraph.Blue,
			asciigraph.Cyan,
			asciigraph.BlueViolet,
			asciigraph.Brown,
			asciigraph.Gray,
			asciigraph.Orange,
			asciigraph.Olive,
		), asciigraph.Height(10)))
	}

	//fmt.Println(asciigraph.PlotMany(p, asciigraph.Precision(0), asciigraph.SeriesColors(
	//	asciigraph.Red,
	//	asciigraph.Yellow,
	//	asciigraph.Green,
	//	asciigraph.Blue,
	//	asciigraph.Cyan,
	//	asciigraph.BlueViolet,
	//	asciigraph.Brown,
	//	asciigraph.Gray,
	//	asciigraph.Orange,
	//	asciigraph.Olive,
	//), asciigraph.Height(10)))
}

func gaussianHist(weight, mean, stdev float64) []float64 {
	y := make([]float64, 100)
	for x := range y {
		y[x] = gaussianY(float64(x), weight, mean, stdev)
	}
	return y
}

func gaussianY(x, a, b, c float64) float64 {
	top := math.Pow(x-b, 2)
	bot := 2 * c * c
	return a * math.Exp(-top/bot)
}

func printLengths(a [][]int) string {
	if len(a) == 0 {
		return ""
	}
	s := new(strings.Builder)
	for i := range a {
		if len(a[i]) == 0 {
			s.WriteString("\tNA")
			continue
		}
		s.WriteString(fmt.Sprintf("\t%d", a[i][0]))
		for j := 1; j < len(a[i]); j++ {
			s.WriteString(fmt.Sprintf(",%d", a[i][j]))
		}
	}
	return s.String()
}

func runMixtureModel(data []int, mm, bestMm *gmm.MixtureModel, f *[]float64) (converged bool, newMm, newBestMm *gmm.MixtureModel) {
	if cap(*f) >= len(data) {
		*f = (*f)[0:len(data)]
	} else {
		*f = make([]float64, len(data))
	}

	for i := range data {
		(*f)[i] = float64(data[i])
	}

	for i := 0; i < 10; i++ {
		converged, _ = gmm.RunMixtureModel(*f, 2, 50, 50, mm)
		if i == 0 {
			mm, bestMm = bestMm, mm
			continue
		}
		if mm.LogLikelihood < bestMm.LogLikelihood {
			mm, bestMm = bestMm, mm
		}
	}
	return converged, mm, bestMm
}

func cleanup(f io.Closer) {
	err := f.Close()
	exception.PanicOnErr(err)
}
