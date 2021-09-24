library(GENOVA)
#load matrix and bed files from HiCPro
WT_10kb <- load_contacts(signal_path = 'AHDC1_WT_10000.matrix', indices_path = 'AHDC1_WT_10000_abs.bed', sample_name = "WT", colour = "blue")
KO_10kb <- load_contacts(signal_path = 'AHDC1_KO_10000.matrix', indices_path = 'AHDC1_KO_10000_abs.bed', sample_name = "KO", colour = "red")
WT_25kb <- load_contacts(signal_path = 'AHDC1_WT_25000.matrix', indices_path = 'AHDC1_WT_25000_abs.bed', sample_name = "WT", colour = "blue")
KO_25kb <- load_contacts(signal_path = 'AHDC1_KO_25000.matrix', indices_path = 'AHDC1_KO_25000_abs.bed', sample_name = "KO", colour = "red")

#load other annotation files of interest
CTCF_down = read.delim('~/Box/Ann_Collier/DiffBind/CTCF_down_hg38.bed', h = F)
CTCF = read.delim('~/Box/Ann_Collier/DiffBind/CTCF_D7_hg38_0.05_IDR.txt.npk.bed', h = F)
GATA3 = read.delim('~/Box/Ann_Collier/GATA3_10k_hg38_0.01.bed', h = F, comment.char = "#")
WT_loops = read.delim('~/Box/Ann_Collier/HiChIP/NovaSeq/AHDC1_WT_10kb.interactions_FitHiC_Q0.01.bed', h= F)
WT_only_loops = read.delim('~/Box/Ann_Collier/HiChIP/NovaSeq/WT_only_10KB_overlap_0.01.bed', h= F)
WT_unique_loops = read.delim('~/Box/Ann_Collier/HiChIP/NovaSeq/WT_unique_loops.bed', h= F)
WT_loops_CTCF_down = read.delim('~/Box/Ann_Collier/HiChIP/AHDC1_MATRIX/WT_loops_CTCF_down.bed', h= F)
WT_loops_GATA3 = read.delim('~/Box/Ann_Collier/HiChIP/AHDC1_MATRIX/WT_loops_GATA3.bed', h= F)
WT_loops_AHDC1_down = read.delim('~/Box/Ann_Collier/HiChIP/AHDC1_MATRIX/WT_loops_AHDC1_down.bed', h= F)
WT_loops_AHDC1_up = read.delim('~/Box/Ann_Collier/HiChIP/AHDC1_MATRIX/WT_loops_AHDC1_up.bed', h= F)
WT_loops_non_AHDC1 = read.delim('~/Box/Ann_Collier/HiChIP/AHDC1_MATRIX/non_AHDC1_target_genes.bed', h= F)

#plot a chromosome matrix
#can overlay with a chip bed file if desired
chromosome_matrix(WT_10kb) #Computes sums of chromosome interactions genome-wide along with expected values, making it easier to identify enriched trans-contacts that may indicate translocations.
hic_matrixplot(exp1 = WT_10kb,exp2 = KO_10kb,
               chrom = 'chr4', start = 100e6, end=112e6, #change coordiantes to location of interest
               loops = list(WT_loops,KO_loops),
               loops.colour = '#998ec3', # purple loops
               loops.type = list('upper','lower'), # only plot in upper triangle
               loops.radius = 20e3, # expand for visibility
               type = 'triangle',
               chip = CTCF_orientation, #optional
               symmAnn = F, # place annotations also on left side
               cut.off = 65) # outer-top

#compute insulation scores
#Insulation scores are calculated by sliding a square along the diagonal of the Hi-C matrix and averaging the values in that square.
insulation = insulation_score(list(WT_10kb, KO_10kb),
                              window = 25)
visualise(insulation, chr = 'chr4', start = 110e6, end=120e6, contrast = 1)
visualise(insulation, contrast = 1)
tornado_insulation(insulation, bed = Lost_insulation, bed_pos = 'center',)
domgram <- insulation_domainogram(list(WT = WT_10kb, KO = KO_10kb),
                                  chrom = "chr4", start = 110e6, end = 120e6) #Creates a domainogram of insulation scores for a genomic region of interest by calculating the insulation scores for a range of sliding square sizes.
#call TADs based on insulation score
#Runs peak detection on the genome wide insulation scores to identify insulation valleys that correspond well to TAD boundaries.
TADcalls <- call_TAD_insulation(insulation)
#overlay TADs with matrix plot
hic_matrixplot(exp1 = WT_10kb,
               chrom = 'chr4',
               start = 110e6,
               end=120e6,
               tads = list(TADcalls$WT, TADcalls$KO), # see ATA
               tads.type = list('lower', 'upper'), # only plot in lower triangle
               tads.colour = c('dodgerblue', 'firebrick2'), # green TAD-borders
               cut.off = 25)
               
#aggregate peak analyis (APA)
#Performs multiple matrix lookup in Hi-C matrices for a twodimensional set of locations, for example loops.
APA_CTCF <- APA(list("WT" = WT_10kb,
                        'KO' = KO_10kb),
                   bedpe = WT_loops_CTCF_down)
visualise(APA_CTCF, metric='diff', colour_lim_contrast = c(-4,4)) #can change the metric to diff or lfc
quantifyAPA_CTCF <- quantify(APA_AHDC1_down, shape="center_vs_quadrants")
boxplot(split(quantifyAPA_CTCF$per_loop$foldchange, f = quantifyAPA_CTCF$per_loop$sample),
        col = c('firebrick2', 'dodgerblue'), outline = F, ylab = 'pixel enrichment at loops')
               
#Get and plot the percentage cis-contacts genome-wide or from a region
cisChrom_out <- cis_trans( list(WT_10kb, KO_10kb) )
barplot(cisChrom_out$cis, names.arg = cisChrom_out$sample, ylim = c(0,100) )
abline(h = 90, col = 'red', lty = 3)
abline(h = 93, col = 'red', lty = 3)

#Relative Contact Probability 
#Produces a dataframe with the probabilities of contacts in a set of distance-bins. Bins are created on a log scale, which leads to equal amounts of datapoints per bin.
RCP_out <- RCP(explist = list(WT_10kb, KO_10kb))
#RCP for specific regions
RCP_beds = RCP(list(WT_10kb, KO_10kb),
              bedlist = list("CTCF" = CTCF,
                             'GATA3' =GATA3,
                             'CTCF_down' =CTCF_down))
visualise(RCP_beds)
visualise(RCP_out)
visualise(RCP_out, smooth = T)
visualise(RCP_beds, contrast = 1, metric = 'lfc')

#compartment score
#need to use at least 25kb resolution
CS_out = compartment_score(list(WT_25kb, KO_25kb), bed = H3K27Ac)
visualise(CS_out, chr = "chr17")
compartment_matrixplot(exp1 = WT_25kb,
                     exp2 = KO_25kb,
                     chrom = 'chr14',
                     arm = 'q',
                     cs.lim = 1.75, # max compartment-score
                     cut.off = 15,
                     chip = H3K27Ac)

#pull out the values between -1 and 1
#plot insulation scores against each other
good_insulation <- ins.score[ins.score$WT > -1 & ins.score$WT <1,]
good_insulation<-na.omit(good_insulation)
plot(x=good_insulation$WT, y=good_insulation$KO, xlim=c(-1, 1), ylim=c(-1,1), ylab="KO insulaiton score", xlab="WT insulation score")
abline(lm(good_insulation$WT~good_insulation$KO), col="red")
ggplot(good_insulation, aes(x=good_insulation$WT, y=good_insulation$KO)) +geom_point(shape=1, alpha=0.1, size = 0.1) +geom_smooth(method=lm, formula = y~x) +ylab("KO inuslation score") +xlab("WT insulation score")  +theme_bw()

#plot comparment scores against each other        
good_compartment <- compartment_scores[compartment_scores$WT > -0.2 & compartment_scores$WT <0.2,]
ggplot(compartment_scores, aes(x=WT, y=KO)) +geom_point(shape=1, alpha=0.1, size = 0.1) +ylab("KO compartment score") +xlab("WT compartment score")  +theme_bw()
ggplot(good_compartment, aes(x=WT, y=KO)) +geom_point(shape=1, alpha=0.1, size = 0.1) +ylab("KO compartment score") +xlab("WT compartment score")  +theme_bw()

#compute enrichment of contacts between TADs
#Calculates the coverage over TADs and between a TAD and its neighbours.
TAD_N_WT <- intra_inter_TAD(list("WT" = WT_10kb,
                                 'KO' = KO_10kb),
                            tad_bed = WT_TADS,
                            max_neighbour = 10)
visualise(TAD_N_WT, geom = 'jitter')
visualise(TAD_N_WT, geom = 'violin', ylim(-5.5))

#aggregate TAD analysis
#Extracts Hi-C matrices around TADs, resizes these to a uniform size and averages the results for all TADs.
ATA <- ATA(list("WT" = WT_10kb,"KO" = KO_10kb), bed = WT_TADS_diff_ins)
visualise(ATA)
visualise(ATA, metric='diff', colour_lim_contrast = c(-5,5), contrast =1)
visualise(ATA_WTcalls,
          colour_lim = c(0,50),
          colour_lim_contrast = c(-5,5),
          focus = 1) 

#directionality index
#The directionality index quantifies the degree of bias between upstream and downstream interactions given a bin on the diagonal. Such biases become apparent near the periphery of TADs: the upstream portion of a TAD interacts more with the downstream bins and inversely, the downstream portion of a TAD interacts more with the upstream bins.
di<-direct_index(list(WT_10kb, KO_10kb), range = 50)
visualise(di, contrast=1, chr = 'chr4', start = 110e6, end=120e6)
di_plot <- di$DI
di_plot<-na.omit(di_plot)
ggplot(di_plot, aes(x=WT, y=KO)) +geom_point(shape=1, alpha=0.1, size = 0.1) +geom_smooth(method=lm, formula = y~x) +ylab("KO directionality index") +xlab("WT directionality index")  +theme_bw()



