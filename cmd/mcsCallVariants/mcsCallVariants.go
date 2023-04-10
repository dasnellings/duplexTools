package main

import (
	"flag"
	"fmt"
	"github.com/dasnellings/MCS_MS/barcode"
	"github.com/dasnellings/MCS_MS/fai"
	"github.com/vertgenlab/gonomics/bed"
	"github.com/vertgenlab/gonomics/cigar"
	"github.com/vertgenlab/gonomics/dna"
	"github.com/vertgenlab/gonomics/exception"
	"github.com/vertgenlab/gonomics/fasta"
	"github.com/vertgenlab/gonomics/fileio"
	"github.com/vertgenlab/gonomics/interval"
	"github.com/vertgenlab/gonomics/sam"
	"github.com/vertgenlab/gonomics/vcf"
	"golang.org/x/exp/slices"
	"io"
	"log"
	"os"
	"runtime"
	"runtime/pprof"
	"sort"
	"strconv"
	"strings"
	"sync"
	"time"
)

func usage() {
	fmt.Print(
		"mcsCallVariants - Call variants from META-CS data processed with annotateReadFamilies.\n" +
			"Usage:\n" +
			"mcsCallVariants [options] -i input.bam -b input.bed -r reference.fasta > output.vcf\n\n")
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
	var excludeBeds inputFiles
	cpuprofile := flag.String("cpuprofile", "", "write cpu profile to `file`")
	memprofile := flag.String("memprofile", "", "write memory profile to `file`")
	input := flag.String("i", "", "Input bam file. Must be indexed.")
	output := flag.String("o", "stdout", "Output VCF file.")
	bedFile := flag.String("b", "", "Input bed file with coordinates of read families, read family ID, and read counts for watson and crick strands. Generated with -bed option in annotateReadFamilies.")
	flag.Var(&excludeBeds, "e", "Bed file(s) with regions to exclude from analysis. May be declared more than once with additional -e flags. Strongly recommended to mask regions with poor mappability.")
	ref := flag.String("r", "", "Fasta file with reference genome used to align input bam. Must be indexed.")
	totalDepth := flag.Int("a", 8, "Minimum total depth of read family for variant consideration.")
	strandedDepth := flag.Int("s", 4, "Minimum depth of independent watson and crick strands for variant consideration")
	endPad := flag.Int("ignoreEnds", 3, "Ignore bases within # of end of a read.")
	minMapQ := flag.Int("minMapQ", 20, "Minimum mapping quality.")
	baseQualPenalty := flag.Float64("baseQualPenalty", 0.25, "Penalty for positions with low quality base. Each reach with a base < minBaseQuality counts towards baseQualPenalty of a read for allele frequency calculations. Note that low quality bases are N-masked and so will always count AGAINST the alternate allele. (e.g. by default each read with a low quality base counts as 0.25 reads for allele frequency determination.")
	allowSuppAln := flag.Bool("allowSupplementaryAlignments", false, "Allow variants using reads that have supplementary alignments annotated.")
	minAf := flag.Float64("minAF", 0.9, "Minimum fraction of reads with alternate allele **Within a read family and within strand** to be considered a variant.")
	minBaseQuality := flag.Int("minBaseQuality", 30, "Minimum base quality to be considered for calling. Bases below threshold will be ignored.")
	maxOverlappingFamilies := flag.Int("maxOverlappingFamilies", 20, "Maximum number of overlapping read families for site to be considered for calling. Low number avoids regions with many misalignments (e.g. centromeres) reducing memory usage. Set to -1 for no limit. Analyzed bed will be `bedfile`.analysis.bed")
	threads := flag.Int("threads", 1, "Number of processor threads to use for calling. Output VCF will be out of order with threads > 1.")
	debugLevel := flag.Int("verbose", 0, "Level of verbosity in log.")
	debugOut := flag.String("debugLog", "", "Print debug logs to file. File may be large. Must be run with threads == 1 for coherent output. ")
	flag.Parse()

	if *cpuprofile != "" {
		f, err := os.Create(*cpuprofile)
		if err != nil {
			log.Fatal("could not create CPU profile: ", err)
		}
		defer f.Close() // error handling omitted for example
		if err := pprof.StartCPUProfile(f); err != nil {
			log.Fatal("could not start CPU profile: ", err)
		}
		defer pprof.StopCPUProfile()
	}

	if len(excludeBeds) == 0 {
		log.Println("WARNING: -e was not declared. It is strongly recommended to mask regions with poor mappability.")
	}

	if *threads == 0 {
		log.Fatal("ERROR: threads must be >= 1.")
	}

	if *input == "" || *bedFile == "" || *ref == "" {
		usage()
		log.Fatal("ERROR: must specify bam (-i), bed (-b), and fasta (-r).")
	}

	if *strandedDepth*2 > *totalDepth {
		log.Fatal("ERROR: -s * 2 should not be larger than -a")
	}

	mcsCallVariants(*input, *output, *ref, *bedFile, excludeBeds, uint8(*minMapQ), *totalDepth, *strandedDepth, *allowSuppAln, *minAf, *minBaseQuality, *baseQualPenalty, *endPad, *maxOverlappingFamilies, *debugLevel, *threads, *debugOut)

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

func mcsCallVariants(input, output, ref, bedFile string, excludeBeds []string, minMapQ uint8, minTotalDepth, minStrandedDepth int, allowSuppAln bool, minAf float64, minBaseQuality int, baseQualPenalty float64, endPad int, maxOverlappingFamilies int, debugLevel int, threads int, debugOut string) {
	// progress tracking
	startTime := time.Now().UnixMilli()

	var excludedRegions map[string]*interval.IntervalNode
	bedFile, excludedRegions = filterInputBed(bedFile, excludeBeds, maxOverlappingFamilies, minTotalDepth, minStrandedDepth)
	vcfOut := fileio.EasyCreate(output)
	vcf.NewWriteHeader(vcfOut, makeVcfHeader(input, ref))
	bedChan := bed.GoReadToChan(bedFile)
	var debugFile io.WriteCloser
	var debugOutChan chan string

	if debugOut != "" {
		debugFile, _ = os.Create(debugOut)
		defer debugFile.Close()
		debugOutChan = make(chan string)
	}

	var err error

	// overhead for multithreading
	wg := new(sync.WaitGroup)
	outputChan := make(chan []vcf.Vcf, 100)
	for i := 0; i < threads; i++ {
		wg.Add(1)
		go spawnThread(bedChan, outputChan, input, ref, minMapQ, minAf, minBaseQuality, baseQualPenalty, endPad, minTotalDepth, minStrandedDepth, allowSuppAln, wg, debugOutChan)
	}

	// spawn a goroutine to wait until threads are done, then close the output
	go func(*sync.WaitGroup) {
		wg.Wait()
		close(outputChan)
		if debugOutChan != nil {
			close(debugOutChan)
		}
	}(wg)

	if debugFile != nil {
		go func() {
			for s := range debugOutChan {
				fmt.Fprintln(debugFile, s)
			}
		}()
	}

	var familiesProcessed int
	var lastVar vcf.Vcf
	lastCheckpointTime := startTime
	currTime := startTime
	for v := range outputChan {
		familiesProcessed++
		if debugLevel > 0 && familiesProcessed%1000 == 0 {
			currTime = time.Now().UnixMilli()
			log.Printf("Processed 1000 Read Families in:\t%dsec\t%s:%d", (currTime-lastCheckpointTime)/1000, lastVar.Chr, lastVar.Pos)
			lastCheckpointTime = currTime
		}
		if len(v) > 0 {
			for i := range v {
				if len(interval.Query(excludedRegions, v[i], "any")) > 0 {
					continue
				}
				vcf.WriteVcf(vcfOut, v[i])
			}
			lastVar = v[len(v)-1]
		}
	}

	endTime := time.Now().UnixMilli()
	log.Printf("Successfully Completed\nRead Families Processed: %d\nTotal Runtime: %d Minutes\n", familiesProcessed, ((endTime-startTime)/1000)/60)

	err = vcfOut.Close()
	exception.PanicOnErr(err)
}

func callFamily(b bed.Bed, bamReader *sam.BamReader, header sam.Header, faSeeker *fasta.Seeker, bai sam.Bai, minMapQ uint8, minAf float64, minBaseQuality int, baseQualPenalty float64, endPad, minTotalDepth, minStrandedDepth int, allowSuppAln bool, recycledReads []sam.Sam, debugOutChan chan<- string) ([]vcf.Vcf, []sam.Sam) {
	var famId string
	var strand byte

	//expectedWatsonDepth, _ := strconv.Atoi(b.Annotation[0])
	//expectedCrickDepth, _ := strconv.Atoi(b.Annotation[1])

	reads := recycledReads[:0]
	reads = sam.SeekBamRegionRecycle(bamReader, bai, b.Chrom, uint32(b.ChromStart), uint32(b.ChromEnd), reads)
	watsonReads := make([]sam.Sam, 0, len(reads))
	crickReads := make([]sam.Sam, 0, len(reads))

	for i := range reads {
		if reads[i].MapQ < minMapQ {
			continue
		}
		sam.ParseExtra(&reads[i])
		famId = barcode.GetRF(&reads[i])
		if famId != b.Name {
			continue
		}
		if hasSuppAln(reads[i]) && !allowSuppAln {
			continue
		}

		clipReadEnds(&reads[i], endPad)
		maskLowQualityBases(&reads[i], minBaseQuality)

		strand = barcode.GetRS(&reads[i])
		if strand == 'W' {
			watsonReads = append(watsonReads, reads[i])
		} else if strand == 'C' {
			crickReads = append(crickReads, reads[i])
		}
	}

	if (len(watsonReads) == 0 && len(crickReads) == 0) || (len(watsonReads) < minStrandedDepth || len(crickReads) < minStrandedDepth) {
		return nil, reads
	}

	sort.Slice(watsonReads, func(i, j int) bool {
		return watsonReads[i].Pos < watsonReads[j].Pos
	})
	sort.Slice(crickReads, func(i, j int) bool {
		return crickReads[i].Pos < crickReads[j].Pos
	})

	watsonPiles := pileup(watsonReads, header)
	crickPiles := pileup(crickReads, header)

	//if debugLevel > 1 && (len(watsonReads) != expectedWatsonDepth || len(crickReads) != expectedCrickDepth) {
	//	log.Printf("WARNING: mismatch in expected (%d/%d) and actual (%d/%d) number of reads, may be supplementary alignments were removed at\n%s\n", expectedWatsonDepth, expectedCrickDepth, len(watsonReads), len(crickReads), b)
	//}

	// remove piles that fall outside the consensus start/end of the read families
	watsonPiles, crickPiles = removePositionalOutliers(watsonPiles, crickPiles, watsonReads, crickReads, endPad, b.Name)
	return pilesToVcfs(watsonPiles, crickPiles, minAf, baseQualPenalty, minStrandedDepth, minTotalDepth, header, faSeeker, b, debugOutChan), reads
}

func pilesToVcfs(watsonPiles, crickPiles []sam.Pile, minAf, baseQualPenalty float64, minStrandedDepth, minTotalDepth int, header sam.Header, faSeeker *fasta.Seeker, b bed.Bed, debugOutChan chan<- string) []vcf.Vcf {
	var variants []vcf.Vcf
	var v vcf.Vcf
	var keeper bool
	var watsonPileIdx, crickPileIdx int

	for { // pos matching between slices of watson and crick piles
		if watsonPileIdx == len(watsonPiles) || crickPileIdx == len(crickPiles) {
			break
		}
		if crickPiles[crickPileIdx].Pos > watsonPiles[watsonPileIdx].Pos {
			watsonPileIdx++
			continue
		}
		if crickPiles[crickPileIdx].Pos < watsonPiles[watsonPileIdx].Pos {
			crickPileIdx++
			continue
		}

		v, keeper = callFromPilePair(watsonPiles[watsonPileIdx], crickPiles[crickPileIdx], minAf, baseQualPenalty, minStrandedDepth, minTotalDepth, header, faSeeker, b, debugOutChan)
		if keeper {
			variants = append(variants, v)
		}

		watsonPileIdx++
		crickPileIdx++
	}
	return variants
}

func callFromPilePair(wPile, cPile sam.Pile, minAf, baseQualPenalty float64, minStrandedDepth, minTotalDepth int, header sam.Header, faSeeker *fasta.Seeker, b bed.Bed, debugOutChan chan<- string) (vcf.Vcf, bool) {
	var watsonDelLen, crickDelLen int
	var watsonInsSeq, crickInsSeq, chr string
	var maxWatsonBase, maxCrickBase dna.Base
	var refBase []dna.Base
	var watsonVarType, crickVarType variantType
	var watsonAltAlleleCount, crickAltAlleleCount, watsonInsAlleleCount, crickInsAlleleCount int
	var err error
	var ans vcf.Vcf

	watsonDepth := pileDepth(wPile, baseQualPenalty)
	crickDepth := pileDepth(cPile, baseQualPenalty)

	//fmt.Printf("evaluating pile %s:%d\nwatson:\t%v\ncrick:\t%v\n\n", header.Chroms[wPile.RefIdx].Name, wPile.Pos, wPile, cPile)

	// matching ref position
	watsonVarType, maxWatsonBase, watsonInsSeq, watsonDelLen, watsonAltAlleleCount, watsonInsAlleleCount = maxBase(wPile)
	crickVarType, maxCrickBase, crickInsSeq, crickDelLen, crickAltAlleleCount, crickInsAlleleCount = maxBase(cPile)

	if debugOutChan != nil {
		debugOutChan <- fmt.Sprintf("watson: %v, crick: %v", wPile, cPile)
	}

	// special case to bias towards insertions since they are assigned to the position before the insertion
	if float64(watsonInsAlleleCount)/float64(watsonDepth) > minAf || float64(crickInsAlleleCount)/float64(crickDepth) > minAf {
		watsonVarType = insertion
		crickVarType = insertion
		watsonAltAlleleCount = watsonInsAlleleCount
		crickAltAlleleCount = crickInsAlleleCount
		if debugOutChan != nil {
			debugOutChan <- fmt.Sprintf("triggered insertion bias")
		}
	}

	// exclude if watson and crick do not agree
	if watsonVarType != crickVarType {
		if debugOutChan != nil {
			debugOutChan <- fmt.Sprintf("variant types do not match, moving on")
		}
		return ans, false
	}

	// exclude if watson or crick AF is less than threshold.
	if float64(watsonAltAlleleCount)/float64(watsonDepth) < minAf || float64(crickAltAlleleCount)/float64(crickDepth) < minAf {
		if debugOutChan != nil {
			debugOutChan <- fmt.Sprintf("does not meet af requirements\nwatson: (%d/%d) = %f\ncrick: (%d/%d) = %f", watsonAltAlleleCount, watsonDepth, float64(watsonAltAlleleCount)/float64(watsonDepth), crickAltAlleleCount, crickDepth, float64(crickAltAlleleCount)/float64(crickDepth))
		}
		return ans, false
	}

	// exclude if below minimum read depth
	if watsonAltAlleleCount < minStrandedDepth || crickAltAlleleCount < minStrandedDepth || watsonAltAlleleCount+crickAltAlleleCount < minTotalDepth {
		if debugOutChan != nil {
			debugOutChan <- fmt.Sprintf("does not meet minimum read depth, moving on")
		}
		return ans, false
	}

	// variant-type specific filters and processing
	chr = header.Chroms[wPile.RefIdx].Name
	switch watsonVarType {
	case snv:
		if maxWatsonBase != maxCrickBase {
			if debugOutChan != nil {
				debugOutChan <- fmt.Sprintf("variant bases do not match, moving on\nwatson: %s\ncrick: %s", dna.BaseToString(maxWatsonBase), dna.BaseToString(maxCrickBase))
			}
			return ans, false
		}

		refBase, err = fasta.SeekByName(faSeeker, chr, int(wPile.Pos-1), int(wPile.Pos))
		dna.AllToUpper(refBase)
		exception.PanicOnErr(err)

		if maxWatsonBase == refBase[0] {
			if debugOutChan != nil {
				debugOutChan <- fmt.Sprintf("alt base matches ref")
			}
			return ans, false
		}
		ans = snvToVcf(wPile, cPile, chr, refBase[0], maxWatsonBase, b.Name)

	case insertion:
		if watsonInsSeq != crickInsSeq {
			if debugOutChan != nil {
				debugOutChan <- fmt.Sprintf("different insertion lengths")
			}
			return ans, false
		}
		ans = insToVcf(wPile, cPile, chr, watsonInsSeq, faSeeker, b.Name)

	case deletion:
		if watsonDelLen != crickDelLen {
			if debugOutChan != nil {
				debugOutChan <- fmt.Sprintf("different deletion lengths")
			}
			return ans, false
		}
		ans = delToVcf(wPile, cPile, chr, watsonDelLen, faSeeker, b.Name)
	}

	return ans, true
}

func pileup(reads []sam.Sam, header sam.Header) []sam.Pile {
	if len(reads) == 0 {
		return nil
	}

	samChan := make(chan sam.Sam, len(reads))
	for i := range reads {
		sclipTerminalIns(&reads[i])
		samChan <- reads[i]
	}
	close(samChan)

	ans := make([]sam.Pile, 0, 100)
	pileChan := sam.GoPileup(samChan, header, false, nil, nil)
	for p := range pileChan {
		ans = append(ans, p)
	}
	return ans
}

type variantType byte

const (
	snv variantType = iota
	insertion
	deletion
	none
)

func maxBase(p sam.Pile) (tp variantType, snvAltBase dna.Base, insSeq string, delLen int, altAlleleCount, maxInsCount int) {
	var maxSnvCount, maxDelCount int

	// check SNV
	for i := 0; i < len(p.CountF); i++ {
		if i == int(dna.Gap) || i == int(dna.N) { // deletions handled below, ignore Ns
			continue
		}
		if p.CountF[i]+p.CountR[i] > maxSnvCount {
			snvAltBase = dna.Base(i)
			maxSnvCount = p.CountF[i] + p.CountR[i]
		}
	}

	// check Del Fwd
	for key := range p.DelCountF {
		if p.DelCountF[key]+p.DelCountR[key] > maxDelCount {
			delLen = key
			maxDelCount = p.DelCountF[key] + p.DelCountR[key]
		}
	}

	// check Del Rev
	for key := range p.DelCountR {
		if p.DelCountF[key]+p.DelCountR[key] > maxDelCount {
			delLen = key
			maxDelCount = p.DelCountF[key] + p.DelCountR[key]
		}
	}

	// check Ins Fwd
	for key := range p.InsCountF {
		if p.InsCountF[key]+p.InsCountR[key] > maxInsCount {
			insSeq = key
			maxInsCount = p.InsCountF[key] + p.InsCountR[key]
		}
	}

	// check Ins Rev
	for key := range p.InsCountR {
		if p.InsCountF[key]+p.InsCountR[key] > maxInsCount {
			insSeq = key
			maxInsCount = p.InsCountF[key] + p.InsCountR[key]
		}
	}

	// score and return winner
	if maxSnvCount > maxInsCount && maxSnvCount > maxDelCount {
		tp = snv
		altAlleleCount = maxSnvCount
		return
	}

	if maxInsCount > maxDelCount {
		tp = insertion
		altAlleleCount = maxInsCount
		return
	}

	if delLen > 0 {
		tp = deletion
		altAlleleCount = maxDelCount
		return
	}

	tp = none
	return
}

func snvToVcf(watsonPile, crickPile sam.Pile, chr string, refBase, altBase dna.Base, readFamily string) vcf.Vcf {
	var v vcf.Vcf
	v.Chr = chr
	v.Pos = int(watsonPile.Pos)
	v.Ref = string(dna.BaseToRune(refBase))
	v.Alt = []string{string(dna.BaseToRune(altBase))}
	v.Filter = "."
	v.Info = "."
	v.Id = "."
	v.Format = []string{"GT", "DP", "WS", "CS", "RF"}

	var totalDepth, watsonDepth, crickDepth string
	totalDepth = fmt.Sprint(calcDepth(watsonPile) + calcDepth(crickPile))
	watsonDepth = fmt.Sprint(watsonPile.CountF[altBase] + watsonPile.CountR[altBase])
	crickDepth = fmt.Sprint(crickPile.CountF[altBase] + crickPile.CountR[altBase])

	v.Samples = make([]vcf.Sample, 1)
	v.Samples[0].Alleles = []int16{1}
	v.Samples[0].FormatData = []string{"", totalDepth, watsonDepth, crickDepth, readFamily}
	return v
}

func insToVcf(watsonPile, crickPile sam.Pile, chr string, insSeq string, faSeeker *fasta.Seeker, readFamily string) vcf.Vcf {
	var v vcf.Vcf
	v.Chr = chr
	v.Pos = int(watsonPile.Pos)

	refBase, err := fasta.SeekByName(faSeeker, chr, int(watsonPile.Pos)-1, int(watsonPile.Pos))
	dna.AllToUpper(refBase)
	exception.PanicOnErr(err)

	v.Ref = string(dna.BaseToRune(refBase[0]))
	v.Alt = []string{string(dna.BaseToRune(refBase[0])) + insSeq}
	v.Filter = "."
	v.Info = "."
	v.Id = "."
	v.Format = []string{"GT", "DP", "WS", "CS", "RF"}

	var totalDepth, watsonDepth, crickDepth string
	totalDepth = fmt.Sprint(calcDepth(watsonPile) + calcDepth(crickPile))
	watsonDepth = fmt.Sprint(watsonPile.InsCountF[insSeq] + watsonPile.InsCountR[insSeq])
	crickDepth = fmt.Sprint(crickPile.InsCountF[insSeq] + crickPile.InsCountR[insSeq])

	v.Samples = make([]vcf.Sample, 1)
	v.Samples[0].Alleles = []int16{1}
	v.Samples[0].FormatData = []string{"", totalDepth, watsonDepth, crickDepth, readFamily}
	return v
}

func delToVcf(watsonPile, crickPile sam.Pile, chr string, delLen int, faSeeker *fasta.Seeker, readFamily string) vcf.Vcf {
	var v vcf.Vcf
	v.Chr = chr
	v.Pos = int(watsonPile.Pos) - 1

	refBase, err := fasta.SeekByName(faSeeker, chr, int(watsonPile.Pos-2), int(watsonPile.Pos-1)+delLen)
	dna.AllToUpper(refBase)
	exception.PanicOnErr(err)

	v.Ref = dna.BasesToString(refBase)
	v.Alt = []string{string(dna.BaseToRune(refBase[0]))}
	v.Filter = "."
	v.Info = "."
	v.Id = "."
	v.Format = []string{"GT", "DP", "WS", "CS", "RF"}

	var totalDepth, watsonDepth, crickDepth string
	totalDepth = fmt.Sprint(calcDepth(watsonPile) + calcDepth(crickPile))
	watsonDepth = fmt.Sprint(watsonPile.DelCountF[delLen] + watsonPile.DelCountR[delLen])
	crickDepth = fmt.Sprint(crickPile.DelCountF[delLen] + crickPile.DelCountR[delLen])

	v.Samples = make([]vcf.Sample, 1)
	v.Samples[0].Alleles = []int16{1}
	v.Samples[0].FormatData = []string{"", totalDepth, watsonDepth, crickDepth, readFamily}
	return v
}

func makeVcfHeader(infile string, referenceFile string) vcf.Header {
	var header vcf.Header
	header.Text = append(header.Text, "##fileformat=VCFv4.2")
	header.Text = append(header.Text, fmt.Sprintf("##reference=%s", referenceFile))
	header.Text = append(header.Text, strings.TrimSuffix(fai.IndexToVcfHeader(fai.ReadIndex(referenceFile+".fai")), "\n"))
	header.Text = append(header.Text, "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">")
	header.Text = append(header.Text, "##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"Total Read Depth\">")
	header.Text = append(header.Text, "##FORMAT=<ID=WS,Number=1,Type=Integer,Description=\"Watson Strand Read Depth\">")
	header.Text = append(header.Text, "##FORMAT=<ID=CS,Number=1,Type=Integer,Description=\"Crick Strand Read Depth\">")
	header.Text = append(header.Text, "##FORMAT=<ID=RF,Number=1,Type=Integer,Description=\"Read Family Identifier\">")
	header.Text = append(header.Text, fmt.Sprintf("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t%s", strings.TrimSuffix(infile, ".bam")))
	return header
}

func removePositionalOutliers(watsonPiles, crickPiles []sam.Pile, watsonReads, crickReads []sam.Sam, endPad int, famId string) (filteredWatsonPiles, filteredCrickPiles []sam.Pile) {
	filteredWatsonPiles = make([]sam.Pile, 0, len(watsonPiles))
	filteredCrickPiles = make([]sam.Pile, 0, len(crickPiles))

	fwdStartMap := make(map[int]int)
	fwdEndMap := make(map[int]int)
	revStartMap := make(map[int]int)
	revEndMap := make(map[int]int)

	for i := range watsonReads {
		if sam.IsPosStrand(watsonReads[i]) {
			fwdStartMap[watsonReads[i].GetChromStart()]++
			fwdEndMap[watsonReads[i].GetChromEnd()]++
		} else {
			revStartMap[watsonReads[i].GetChromStart()]++
			revEndMap[watsonReads[i].GetChromEnd()]++
		}
	}
	for i := range crickReads {
		if sam.IsPosStrand(crickReads[i]) {
			fwdStartMap[crickReads[i].GetChromStart()]++
			fwdEndMap[crickReads[i].GetChromEnd()]++
		} else {
			revStartMap[crickReads[i].GetChromStart()]++
			revEndMap[crickReads[i].GetChromEnd()]++
		}
	}

	var fwdStart, fwdEnd, revStart, revEnd, maxCount int
	for key, val := range fwdStartMap {
		if val > maxCount || (val == maxCount && key < fwdStart) {
			fwdStart = key
			maxCount = val
		}
	}
	maxCount = 0
	for key, val := range fwdEndMap {
		if val > maxCount || (val == maxCount && key > fwdEnd) {
			fwdEnd = key
			maxCount = val
		}
	}
	maxCount = 0
	for key, val := range revStartMap {
		if val > maxCount || (val == maxCount && key < revStart) {
			revStart = key
			maxCount = val
		}
	}
	maxCount = 0
	for key, val := range revEndMap {
		if val > maxCount || (val == maxCount && key > revEnd) {
			revEnd = key
			maxCount = val
		}
	}

	for i := range watsonPiles {
		if (int(watsonPiles[i].Pos) > fwdStart && int(watsonPiles[i].Pos) < fwdEnd) ||
			(int(watsonPiles[i].Pos) > revStart && int(watsonPiles[i].Pos) < revEnd) {
			filteredWatsonPiles = append(filteredWatsonPiles, watsonPiles[i])
		}
	}

	for i := range crickPiles {
		if (int(crickPiles[i].Pos) > fwdStart && int(crickPiles[i].Pos) < fwdEnd) ||
			(int(crickPiles[i].Pos) > revStart && int(crickPiles[i].Pos) < revEnd) {
			filteredCrickPiles = append(filteredCrickPiles, crickPiles[i])
		}
	}
	return
}

// calcDepth returns the number of reads in the input pile
func calcDepth(s sam.Pile) int {
	var depth int
	for i := range s.CountF {
		depth += s.CountF[i] + s.CountR[i]
	}
	return depth
}

// sclipTerminalIns will convert an insertion on the left or right end of the read to a soft clip
func sclipTerminalIns(s *sam.Sam) {
	if len(s.Cigar) == 0 || s.Cigar[0].Op == '*' {
		return
	}
	if s.Cigar[0].Op == 'I' {
		s.Cigar[0].Op = 'S'
	}
	if s.Cigar[len(s.Cigar)-1].Op == 'I' {
		s.Cigar[len(s.Cigar)-1].Op = 'S'
	}

	// catch case where beginning/end of read is already soft clipped
	if len(s.Cigar) >= 2 && s.Cigar[0].Op == 'S' && s.Cigar[1].Op == 'I' {
		s.Cigar[1].Op = 'S'
		s.Cigar[1].RunLength += s.Cigar[0].RunLength
		s.Cigar = s.Cigar[1:]
	}

	if len(s.Cigar) >= 2 && s.Cigar[len(s.Cigar)-1].Op == 'S' && s.Cigar[len(s.Cigar)-2].Op == 'I' {
		s.Cigar[len(s.Cigar)-2].Op = 'S'
		s.Cigar[len(s.Cigar)-2].RunLength += s.Cigar[len(s.Cigar)-1].RunLength
		s.Cigar = s.Cigar[:len(s.Cigar)-1]
	}
}

func filterInputBed(bedFile string, excludeBeds []string, maxOverlaps int, minTotalDepth int, minStrandedDepth int) (string, map[string]*interval.IntervalNode) {
	var excludeIntervals []interval.Interval
	var tree map[string]*interval.IntervalNode
	for _, e := range excludeBeds {
		bChan := bed.GoReadToChan(e)
		for b := range bChan {
			excludeIntervals = append(excludeIntervals, b)
		}
	}
	tree = interval.BuildTree(excludeIntervals)

	outfile := strings.TrimSuffix(bedFile, ".bed") + ".analysis.bed"
	beds := bed.GoReadToChan(bedFile)
	out := fileio.EasyCreate(outfile)
	overlaps := make([]bed.Bed, 0, 1000)
	var watsonDepth, crickDepth int
	for b := range beds {
		switch {
		case len(overlaps) == 0:
			overlaps = append(overlaps, b)

		case bed.Overlap(overlaps[0], b):
			overlaps = append(overlaps, b)

		default: // does not overlap
			if len(overlaps) <= maxOverlaps { // write
				for i := range overlaps {
					watsonDepth, _ = strconv.Atoi(overlaps[i].Annotation[0])
					crickDepth, _ = strconv.Atoi(overlaps[i].Annotation[1])
					if watsonDepth+crickDepth < minTotalDepth {
						continue
					}
					if watsonDepth < minStrandedDepth || crickDepth < minStrandedDepth {
						continue
					}
					if len(excludeBeds) > 0 && len(interval.Query(tree, overlaps[i], "di")) > 0 { // query entirely contained within excluded region
						continue
					}
					bed.WriteBed(out, overlaps[i])
				}
			}
			overlaps = overlaps[:0]
			overlaps = append(overlaps, b)
		}
	}

	if len(overlaps) == 1 {
		bed.WriteBed(out, overlaps[0])
	}
	err := out.Close()
	exception.PanicOnErr(err)
	return outfile, tree
}

func spawnThread(inputChan <-chan bed.Bed, outputChan chan<- []vcf.Vcf, inputBam string, ref string, minMapQ uint8, minAf float64, minBaseQuality int, baseQualPenalty float64, endPad, minTotalDepth, minStrandedDepth int, allowSuppAln bool, wg *sync.WaitGroup, debugOutChan chan<- string) {
	bamReader, bamHeader := sam.OpenBam(inputBam)
	bai := sam.ReadBai(inputBam + ".bai")
	faSeeker := fasta.NewSeeker(ref, "")
	var err error

	var familyVariants []vcf.Vcf
	var recycledReads []sam.Sam
	for b := range inputChan {
		familyVariants, recycledReads = callFamily(b, bamReader, bamHeader, faSeeker, bai, minMapQ, minAf, minBaseQuality, baseQualPenalty, endPad, minTotalDepth, minStrandedDepth, allowSuppAln, recycledReads, debugOutChan)
		outputChan <- familyVariants
	}

	err = bamReader.Close()
	exception.PanicOnErr(err)
	err = faSeeker.Close()
	exception.PanicOnErr(err)
	wg.Done()
}

func clipReadEnds(s *sam.Sam, clipLen int) {
	if s.Cigar == nil || len(s.Cigar) == 0 || s.Cigar[0].Op == '*' {
		return
	}

	var anyNonClip bool
	for i := range s.Cigar {
		if s.Cigar[i].Op != 'S' {
			anyNonClip = true
			break
		}
	}

	if !anyNonClip {
		return
	}

	clipFwd(s, clipLen)
	clipRev(s, clipLen)

	// collapse cigar if everything is soft clipped
	if len(s.Cigar) == 2 && s.Cigar[0].Op == 'S' && s.Cigar[1].Op == 'S' {
		s.Cigar[0].RunLength += s.Cigar[1].RunLength
		s.Cigar = s.Cigar[:1]
	}

	//if cigar.QueryLength(s.Cigar) != len(s.Seq) {
	//	log.Panic("something went horribly wrong with cigar\n", s)
	//}
}

func clipFwd(s *sam.Sam, clipLen int) {
	if clipLen < 1 {
		return
	}

	// check if first index is soft clip, if not make a soft clip with len = 0
	if s.Cigar[0].Op != 'S' {
		s.Cigar = slices.Insert(s.Cigar, 0, cigar.Cigar{Op: 'S', RunLength: 0})
	}
	var numToClip int = clipLen
	var currNumToClip int
	for i := 1; numToClip > 0; i++ {
		// increment pos as well as cigar
		switch s.Cigar[i].Op {
		case 'M':
			currNumToClip = min(s.Cigar[i].RunLength, numToClip)
			s.Cigar[i].RunLength -= currNumToClip
			s.Cigar[0].RunLength += currNumToClip
			s.Pos += uint32(currNumToClip)
			numToClip -= currNumToClip

		case 'D':
			s.Pos += uint32(s.Cigar[i].RunLength)
			s.Cigar[i].RunLength = 0

		case 'I':
			currNumToClip = min(s.Cigar[i].RunLength, numToClip)
			s.Cigar[0].RunLength += currNumToClip
			s.Cigar[i].RunLength -= currNumToClip
			numToClip -= currNumToClip

		case 'S':
			s.Cigar = cleanCigar(s.Cigar)
			return
		}
	}
	s.Cigar = cleanCigar(s.Cigar)
}

func clipRev(s *sam.Sam, clipLen int) {
	if clipLen < 1 {
		return
	}

	// check if last index is soft clip, if not make a soft clip with len = 0
	if s.Cigar[len(s.Cigar)-1].Op != 'S' {
		s.Cigar = append(s.Cigar, cigar.Cigar{Op: 'S', RunLength: 0})
	}
	var numToClip int = clipLen
	var currNumToClip int
	lastIdx := len(s.Cigar) - 1
	for i := lastIdx - 1; numToClip > 0; i-- {
		// increment pos as well as cigar
		switch s.Cigar[i].Op {
		case 'M', 'I':
			currNumToClip = min(s.Cigar[i].RunLength, numToClip)
			s.Cigar[i].RunLength -= currNumToClip
			s.Cigar[lastIdx].RunLength += currNumToClip
			numToClip -= currNumToClip

		case 'D':
			s.Cigar[i].RunLength = 0

		case 'S':
			s.Cigar = cleanCigar(s.Cigar)
			return
		}
	}
	s.Cigar = cleanCigar(s.Cigar)
}

func pileDepth(p sam.Pile, baseQualPenalty float64) int {
	var depth int
	var maskCount int
	for i := range p.CountF {
		if i == int(dna.N) {
			maskCount++
			continue
		}
		depth += p.CountF[i] + p.CountR[i]
	}
	depth += int(float64(maskCount) * baseQualPenalty)
	return depth
}

func maskLowQualityBases(s *sam.Sam, minQual int) {
	var currQual uint8
	for i := range s.Qual {
		currQual = s.Qual[i] - 33
		if currQual < uint8(minQual) {
			s.Seq[i] = dna.N
		}
	}
}

func cleanCigar(c []cigar.Cigar) []cigar.Cigar {
	// remove all indexes with RunLength of 0
	for i := 0; i < len(c); i++ {
		if c[i].RunLength == 0 {
			c = slices.Delete(c, i, i+1)
			i--
		}
	}
	return c
}

func hasSuppAln(r sam.Sam) bool {
	_, found, err := sam.QueryTag(r, "SA")
	if err != nil || !found {
		return false
	}
	return true
}

func min(a, b int) int {
	if a < b {
		return a
	}
	return b
}

func max(a, b int) int {
	if a > b {
		return a
	}
	return b
}
