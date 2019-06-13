# NOTE: HOMER PATHS DON'T WORK ANYMORE DUE TO ORGANIZATION 6/11/19
# ==========================================================
# QUESTION 2_1
# ==========================================================
# What % of eQTLs (cis / trans consider separately) overlap 
# a motif or are in LD with a SNP that overlaps a motif? 
# If a SNP is not genotyped or poorly imputed (missing data or incorrect data), 
# we won’t be able to see it’s true association with gene expression. 
# Therefore, a SNP in tight LD with this poorly documented SNP (which is in reality the causal SNP) 
# might show an artificially strong association signal with expression. 

# # use homer to find trans eQTLs that also overlap cis eQTLs
# write.table(franke_cis_trans_common_snps, file="q2_1_written/franke_cis_trans_common.txt",row.names= FALSE,col.names=FALSE,quote=FALSE) # write SNPs to file
# 
# # match cis/trans snps to cis genes, keep only unique genes for Homer findmotifs.pl // do I need more specificity for gene?  i.e. gene.1, gene.2...?  Do i keep only unique genes?
# franke_cis_trans_common_cis_genes <- unique(data.table(right_join(franke_cis_gene_snps, franke_cis_trans_common_snps))[,"Gene"]) # 24490 merged -> 4605 unique
# write.table(franke_cis_trans_common_cis_genes, file="q2_1_input/homer/franke_cis_trans_common_cis_genes.txt",row.names= FALSE,col.names="Gene_ID",quote=FALSE) # write cis genes to file
# # match cis/trans snps to trans genes, keep only unique genes for Homer findmotifs.pl // same questions as for cis
# franke_cis_trans_common_trans_genes <- unique(data.table(right_join(data.table(franke_trans_data[,"SNP"], franke_trans_data[,"Gene"]), franke_cis_trans_common_snps))[,"Gene"]) #56038 merged -> 6059 unique
# write.table(franke_cis_trans_common_trans_genes, file="q2_1_input/homer/franke_cis_trans_common_trans_genes.txt",row.names= FALSE,col.names="Gene_ID",quote=FALSE) # write trans genes to file
# 
# # having issues getting homer to find file?  am i putting it in the right directory? FIXED - spell "franke" correctly
# # findMotifs.pl franke_cis_trans_common_trans_genes.txt human  motifResults_trans1/ -find data/knownTFs/vertebrates/known.motifs > /my_dir/trans_output.txt
# 
# trans_homer_results <- fread("q2_1_input/homer/trans_output.txt",header = TRUE, sep = "\t", dec = ".") #336,445
#   
# # parse out after parentheses in motif name
# # create table from data + trans_genes
# # print list of factors
# 
# trans_homer_results$MotifNameAbbreviated <-  sub("*\\(.*", "", trans_homer_results$'Motif Name')
# # trans_homer_results_organized <- data.table(trans_homer_results[,"Offset"],trans_homer_results[,"Strand"],trans_homer_results[,"MotifScore"],trans_homer_results[,"Ensembl"],trans_homer_results[,"MotifNameAbbreviated"]) #336,445
# # franke_trans_gene_snps <- data.table(franke_trans_data$SNP, franke_trans_data$Gene) #59,786
# # colnames(franke_trans_gene_snps) <- c("SNP", "Ensembl")
# # trans_homer_results_organized_snps <- left_join(trans_homer_results_organized, franke_trans_gene_snps) #left_join only yields 3,359,095, inner_join only yields 3,359,095 results
# # colnames(trans_homer_results_organized_snps) <- c("Offset", "Strand", "MotifScore", "Trans_Gene", "MotifNameAbbreviated", "SNP")
# # trans_homer_results_organized_snps_cis_gene_snps <- left_join(trans_homer_results_organized_snps, franke_cis_gene_snps) #left_join only yields 30,472,391, inner_join only yields 30,118,485
# # colnames(trans_homer_results_organized_snps_cis_gene_snps) <- c("Offset", "Strand", "MotifScore", "Trans_Gene", "MotifNameAbbreviated", "SNP", "Cis_Gene")
# 
# # WITHOUT STRAND
# trans_homer_results_organized <- data.table(trans_homer_results[,"Offset"],trans_homer_results[,"MotifScore"],trans_homer_results[,"Ensembl"],trans_homer_results[,"MotifNameAbbreviated"]) #336,445
# franke_trans_gene_snps <- data.table(franke_trans_data$SNP, franke_trans_data$Gene) #59,786
# colnames(franke_trans_gene_snps) <- c("SNP", "Ensembl")
# trans_homer_results_organized_snps <- left_join(trans_homer_results_organized, franke_trans_gene_snps) #left_join only yields 3,359,095, inner_join only yields 3,359,095 results
# colnames(trans_homer_results_organized_snps) <- c("Offset", "MotifScore", "Trans_Gene", "MotifNameAbbreviated", "SNP")
# trans_homer_results_organized_snps_cis_gene_snps <- left_join(trans_homer_results_organized_snps, franke_cis_gene_snps) #left_join only yields 30,472,391, inner_join only yields 30,118,485
# colnames(trans_homer_results_organized_snps_cis_gene_snps) <- c("Offset", "MotifScore", "Trans_Gene", "MotifNameAbbreviated", "SNP", "Cis_Gene")
# 
# write.csv(trans_homer_results_organized_snps_cis_gene_snps, file="q2_1_input/homer/trans_homer_results_organized_snps_cis_gene_snps.csv")
# # WRITE CODE TO GET UNIQUE ONES - can't run on local computer
# 
# franke_cis_trans_common_snps_gene <- inner_join(franke_trans_gene_snps, franke_cis_gene_snps, by="SNP")

# ==========================================================
# QUESTION 2-1: TIFFANY'S SCRIPT FOR INTERESTING SNPS
# ==========================================================
# 
# franke_cis_gene_snps_zscore <- data.table(franke_cis_data$SNP, franke_cis_data$Gene, franke_cis_data$Zscore) 
# colnames(franke_cis_gene_snps_zscore) <- c("SNP", "Ensembl", "Zscore")
# franke_trans_gene_snps_zscore <- data.table(franke_trans_data$SNP, franke_trans_data$Gene, franke_trans_data$Zscore)
# colnames(franke_trans_gene_snps_zscore) <- c("SNP", "Ensembl", "Zscore")
# franke_cis_trans_common_snps_gene_zscore <- inner_join(franke_cis_gene_snps_zscore, franke_trans_gene_snps_zscore, by="SNP")
# 
# #setup: restricted ourselves to SNPs that are both cis and trans eQTLs.
# #did not restrict ourselves to cis Genes that are TFs.
# #goal: want to identify if any cis eQTL Genes are TFs that can bind in the promoter of the trans eQTL gene
# dat <- fread("q2_1_input/homer/trans_output.txt", header = T, sep = "\t", dec = ".")
# #336K by 16
# 
# #table with SNP, cis Gene, trans Gene, list of Motifs
# TFs <- sapply(1:nrow(dat), function(x) strsplit(dat$`Motif Name`[x], split = "[(]")[[1]][1])
# length(unique(TFs)) #407
# 
# #snps <- fread("q2_1_input/franke_cis_trans_common_snps_gene.csv", header = T, sep = ",")
# snps <- franke_cis_trans_common_snps_gene_zscore
# colnames(snps) <- c("SNP","transGene","transZscore", "cisGene", "cisZscore")
# # 601288  x    5
# 
# # franke_cis_data <- fread("input_data/cis-eQTL_significant_20181017.txt", header = T, sep = "\t", dec = ".")
# common_cisgenename <- franke_cis_data$GeneSymbol[match(snps$cisGene, franke_cis_data$Gene)]
# 
# #did not restrict ourselves to cis Genes that are TFs.
# #the only possibility for explaining this first regime is if the cis eQTL gene is a TF gene
# m <- match(toupper(common_cisgenename), toupper(unique(TFs)))
# #important to match case! # UPDATED LIST: these numbers don't match
# length(which(is.na(m))) #590719, NOW 590847
# length(which(!is.na(m))) #10569 (can reduce our search space to this), NOW 10441
# 
# w <- which(!is.na(m))
# snps_interesting <- snps[w,]
# common_cisgenename_interesting <- common_cisgenename[w]
# #sanity check
# z <- match(toupper(common_cisgenename_interesting), toupper(unique(TFs)))
# length(which(is.na(z))) == 0
# 
# #########
# #motif_matches <- sapply(1:nrow(snps), function(x) paste(TFs[which(dat$Ensembl == as.character(snps[x,3]))], collapse = ","))
# matches <- c()
# for (i in 1:nrow(snps_interesting)){
#   if (i %% 1000 == 0){print(i)}
#   matches[i] <- match(toupper(common_cisgenename_interesting[i]), toupper(TFs[which(dat$Ensembl == as.character(snps_interesting[i,"transGene"]))]))
# }
# length(which(!is.na(matches))) #458, NOW 812
# 
# #snps_interesting_TFmatch <- cbind(snps_interesting, ifelse(is.na(matches), 0, 1))
# #w <- which(snps_interesting_TFmatch$V2 == 1)
# snps_interesting_TFmatch <- snps_interesting[which(!is.na(matches)),]
# snps_interesting_TFmatch <- cbind(snps_interesting_TFmatch, common_cisgenename_interesting[which(!is.na(matches))])
# colnames(snps_interesting_TFmatch)[6] <- "cisGene_commonname"
# #snps_interesting_TFmatch <- snps_interesting_TFmatch[,c(1,2,4,3)]
# write.table(snps_interesting_TFmatch, "q2_1_written/snps_interesting_TFmatch.txt", row.names = F, col.names = T)
# 
# ==========================================================
# ACTIVATOR/REPRESSOR ANALYSIS AND COMPARISON WITH Z-SCORES
# ==========================================================

snps_interesting_TFmatch <- fread("q2_1_written/snps_interesting_TFmatch.txt")
# A represents Activator, R represents Repressor, A/R indicates it can be both.  Information taken from UniProt
z <- snps_interesting_TFmatch[!duplicated(snps_interesting_TFmatch$"cisGene_commonname"),2:6]
cis_gene_cards <- integer(28)
for (i in 1:nrow(z)){
  cis_gene_cards[i] <- paste("https://www.genecards.org/cgi-bin/carddisp.pl?gene=", z[i,"cisGene_commonname"], sep="")
}

functions_u <- data.frame(
  "A - Transcriptional activator (By similarity). Binds a GC box motif. Could play a role in B-cell growth and development.", 
  "A/R - Nuclear receptor that binds DNA as a monomer to ROR response elements (RORE) containing a single core motif half-site 5'-AGGTCA-3' preceded by a short A-T-rich sequence. Key regulator of embryonic development, cellular differentiation, immunity, circadian rhythm as well as lipid, steroid, xenobiotics and glucose metabolism. Considered to have intrinsic transcriptional activity, have some natural ligands like oxysterols that act as agonists (25-hydroxycholesterol) or inverse agonists (7-oxygenated sterols), enhancing or repressing the transcriptional activity, respectively. Recruits distinct combinations of cofactors to target genes regulatory regions to modulate their transcriptional expression, depending on the tissue, time and promoter contexts. Regulates genes involved in photoreceptor development including OPN1SW, OPN1SM and ARR3 and skeletal muscle development with MYOD1. Required for proper cerebellum development (PubMed:29656859). Regulates SHH gene expression, among others, to induce granule cells proliferation as well as expression of genes involved in calcium-mediated signal transduction. Regulates the circadian expression of several clock genes, including CLOCK, ARNTL/BMAL1, NPAS2 and CRY1. Competes with NR1D1 for binding to their shared DNA response element on some clock genes such as ARNTL/BMAL1, CRY1 and NR1D1 itself, resulting in NR1D1-mediated repression or RORA-mediated activation of clock genes expression, leading to the circadian pattern of clock genes expression. Therefore influences the period length and stability of the clock. Regulates genes involved in lipid metabolism such as apolipoproteins APOA1, APOA5, APOC3 and PPARG. In liver, has specific and redundant functions with RORC as positive or negative modulator of expression of genes encoding phase I and phase II proteins involved in the metabolism of lipids, steroids and xenobiotics, such as CYP7B1 and SULT2A1. Induces a rhythmic expression of some of these genes. In addition, interplays functionally with NR1H2 and NR1H3 for the regulation of genes involved in cholesterol metabolism. Also involved in the regulation of hepatic glucose metabolism through the modulation of G6PC and PCK1. In adipose tissue, plays a role as negative regulator of adipocyte differentiation, probably acting through dual mechanisms. May suppress CEBPB-dependent adipogenesis through direct interaction and PPARG-dependent adipogenesis through competition for DNA-binding. Downstream of IL6 and TGFB and synergistically with RORC isoform 2, is implicated in the lineage specification of uncommitted CD4(+) T-helper (T(H)) cells into T(H)17 cells, antagonizing the T(H)1 program. Probably regulates IL17 and IL17F expression on T(H) by binding to the essential enhancer conserved non-coding sequence 2 (CNS2) in the IL17-IL17F locus. Involved in hypoxia signaling by interacting with and activating the transcriptional activity of HIF1A. May inhibit cell growth in response to cellular stress. May exert an anti-inflammatory role by inducing CHUK expression and inhibiting NF-kappa-B signaling. RORA_HUMAN,P35398", 
  "A - Transcription activator which binds specifically to the MEF2 element present in the regulatory regions of many muscle-specific genes. Controls cardiac morphogenesis and myogenesis, and is also involved in vascular development. Plays an essential role in hippocampal-dependent learning and memory by suppressing the number of excitatory synapses and thus regulating basal and evoked synaptic transmission. Crucial for normal neuronal development, distribution, and electrical activity in the neocortex. Necessary for proper development of megakaryocytes and platelets and for bone marrow B-lymphopoiesis. Required for B-cell survival and proliferation in response to BCR stimulation, efficient IgG1 antibody responses to T-cell-dependent antigens and for normal induction of germinal center B-cells. May also be involved in neurogenesis and in the development of cortical architecture (By similarity). Isoform 3 and isoform 4, which lack the repressor domain, are more active than isoform 1 and isoform 2.", 
  "Unknown - Binds to the CACCC box of erythroid cell-expressed genes. May play a role in hematopoiesis", 
  "A - Signal transducer and transcription activator that mediates cellular responses to interferons (IFNs), cytokine KITLG/SCF and other cytokines and other growth factors. Following type I IFN (IFN-alpha and IFN-beta) binding to cell surface receptors, signaling via protein kinases leads to activation of Jak kinases (TYK2 and JAK1) and to tyrosine phosphorylation of STAT1 and STAT2. The phosphorylated STATs dimerize and associate with ISGF3G/IRF-9 to form a complex termed ISGF3 transcription factor, that enters the nucleus (PubMed:28753426). ISGF3 binds to the IFN stimulated response element (ISRE) to activate the transcription of IFN-stimulated genes (ISG), which drive the cell in an antiviral state. In response to type II IFN (IFN-gamma), STAT1 is tyrosine- and serine-phosphorylated (PubMed:26479788). It then forms a homodimer termed IFN-gamma-activated factor (GAF), migrates into the nucleus and binds to the IFN gamma activated sequence (GAS) to drive the expression of the target genes, inducing a cellular antiviral state. Becomes activated in response to KITLG/SCF and KIT signaling. May mediate cellular responses to activated FGFR1, FGFR2, FGFR3 and FGFR4. STAT1_HUMAN,P42224", 
  "A - Transcriptional activator; DNA-binding protein that specifically recognize the sequence 5'-YAAC[GT]G-3'. Plays an important role in the control of proliferation and differentiation of hematopoietic progenitor cells. MYB_HUMAN,P10242", 
  "R - Transcription factor associated with the BAF SWI/SNF chromatin remodeling complex (By similarity). Repressor of fetal hemoglobin (HbF) level (PubMed:26375765). Involved in brain development (PubMed:27453576). Functions as a myeloid and B-cell proto-oncogene. May play important roles in leukemogenesis and hematopoiesis. Essential factor in lymphopoiesis required for B-cell formation in fetal liver. May function as a modulator of the transcriptional repression activity of ARP1 (By similarity).", 
  "A - Sequence-specific DNA-binding transcription factor.  Binds to two specific DNA sites located in the promoter region of HOXA4.", 
  "A/R - Transcription regulator. Forms a sequence-specific DNA-binding protein complex with MYC or MAD which recognizes the core sequence 5'-CAC[GA]TG-3'. The MYC:MAX complex is a transcriptional activator, whereas the MAD:MAX complex is a repressor. May repress transcription via the recruitment of a chromatin remodeling complex containing H3 'Lys-9' histone methyltransferase activity. Represses MYC transcriptional activity from E-box elements.", 
  "R - Receptor for retinoic acid. Retinoic acid receptors bind as heterodimers to their target response elements in response to their ligands, all-trans or 9-cis retinoic acid, and regulate gene expression in various biological processes. The RXR/RAR heterodimers bind to the retinoic acid response elements (RARE) composed of tandem 5'-AGGTCA-3' sites known as DR1-DR5. In the absence of ligand, the RXR-RAR heterodimers associate with a multiprotein complex containing transcription corepressors that induce histone acetylation, chromatin condensation and transcriptional suppression. On ligand binding, the corepressors dissociate from the receptors and associate with the coactivators leading to transcriptional activation. RARA plays an essential role in the regulation of retinoic acid-induced germ cell development during spermatogenesis. Has a role in the survival of early spermatocytes at the beginning prophase of meiosis. In Sertoli cells, may promote the survival and development of early meiotic prophase spermatocytes. In concert with RARG, required for skeletal growth, matrix homeostasis and growth plate function (By similarity). RARA_HUMAN,P10276", 
  "R - Transcriptional repressor involved in the regulation of the circadian rhythm by negatively regulating the activity of the clock genes and clock-controlled genes. Acts as the negative limb of a novel autoregulatory feedback loop (DEC loop) which differs from the one formed by the PER and CRY transcriptional repressors (PER/CRY loop). Both these loops are interlocked as it represses the expression of PER1 and in turn is repressed by PER1/2 and CRY1/2. Represses the activity of the circadian transcriptional activator: CLOCK-ARNTL/BMAL1 heterodimer by competing for the binding to E-box elements (5'-CACGTG-3') found within the promoters of its target genes. Negatively regulates its own expression and the expression of DBP and BHLHE41/DEC2. Acts as a corepressor of RXR and the RXR-LXR heterodimers and represses the ligand-induced RXRA/B/G, NR1H3/LXRA, NR1H4 and VDR transactivation activity. Inhibits HNF1A-mediated transactivation of CYP1A2, CYP2E1 AND CYP3A11 (By similarity).", 
  "A/R - Transcription factor; can act both as activator and as repressor. Binds the 5'-CACCC-3' core sequence. Binds to the promoter region of its own gene and can activate its own transcription. Regulates the expression of key transcription factors during embryonic development. Plays an important role in maintaining embryonic stem cells, and in preventing their differentiation. Required for establishing the barrier function of the skin and for postnatal maturation and maintenance of the ocular surface. Involved in the differentiation of epithelial cells and may also function in skeletal and kidney development. Contributes to the down-regulation of p53/TP53 transcription.", 
  "A - Transcriptional activator. Binds to the interferon-stimulated response element (ISRE) of the MHC class I promoter. Binds the immunoglobulin lambda light chain enhancer, together with PU.1. Probably plays a role in ISRE-targeted signal transduction mechanisms specific to lymphoid cells. Involved in CD8(+) dendritic cell differentiation by forming a complex  with the BATF-JUNB heterodimer in immune cells, leading to recognition of AICE sequence (5'-TGAnTCA/GAAA-3'),  an immune-specific regulatory element, followed by cooperative binding of BATF and IRF4 and activation of genes  (By similarity)", 
  "A - Acts as a transcriptional regulator of PAX6. Acts as a transcriptional activator of PF4 in complex with PBX1 or PBX2. Required for hematopoiesis, megakaryocyte lineage development and vascular patterning. May function as a cofactor for HOXA7 and HOXA9 in the induction of myeloid leukemias. MEIS1_HUMAN,O00470",
  "R - Acts as a transcriptional repressor. Inhibits interleukin-2 (IL-2) gene expression. Enhances or represses the promoter activity of the ATP1A1 gene depending on the quantity of cDNA and on the cell type. Represses E-cadherin promoter and induces an epithelial-mesenchymal transition (EMT) by recruiting SMARCA4/BRG1. Represses BCL6 transcription in the presence of the corepressor CTBP1. Positively regulates neuronal differentiation. Represses RCOR1 transcription activation during neurogenesis. Represses transcription by binding to the E box (5'-CANNTG-3'). Promotes tumorigenicity by repressing stemness-inhibiting microRNAs. ZEB1_HUMAN,P37275",
  "R - Transcription factor that is the main target of insulin signaling and regulates metabolic homeostasis in response to oxidative stress. Binds to the insulin response element (IRE) with consensus sequence 5'-TT[G/A]TTTTG-3' and the related Daf-16 family binding element (DBE) with consensus sequence 5'-TT[G/A]TTTAC-3'. Activity suppressed by insulin. Main regulator of redox balance and osteoblast numbers and controls bone mass. Orchestrates the endocrine function of the skeleton in regulating glucose metabolism. Acts synergistically with ATF4 to suppress osteocalcin/BGLAP activity, increasing glucose levels and triggering glucose intolerance and insulin insensitivity. Also suppresses the transcriptional activity of RUNX2, an upstream activator of osteocalcin/BGLAP. In hepatocytes, promotes gluconeogenesis by acting together with PPARGC1A and CEBPA to activate the expression of genes such as IGFBP1, G6PC and PCK1. Important regulator of cell death acting downstream of CDK1, PKB/AKT1 and STK4/MST1. Promotes neural cell death. Mediates insulin action on adipose tissue. Regulates the expression of adipogenic genes such as PPARG during preadipocyte differentiation and, adipocyte size and adipose tissue-specific gene expression in response to excessive calorie intake. Regulates the transcriptional activity of GADD45A and repair of nitric oxide-damaged DNA in beta-cells. Required for the autophagic cell death induction in response to starvation or oxidative stress in a transcription-independent manner. Mediates the function of MLIP in cardiomyocytes hypertrophy and cardiac remodeling (By similarity).",
  "A - Functions as a transcriptional activator playing a crucial role during development. Functions in trophoblast differentiation and later in gastrulation, regulating both mesoderm delamination and endoderm specification. Plays a role in brain development being required for the specification and the proliferation of the intermediate progenitor cells and their progeny in the cerebral cortex. Also involved in the differentiation of CD8+ T-cells during immune response regulating the expression of lytic effector genes. EOMES_HUMAN,O95936", 
  "A - Transcriptional activator which recognizes variations of the palindromic sequence 5'-ATTCCCNNGGGAATT-3'.",
  "R - Acts as a transcriptional regulator that recognizes and binds to the sequence 5'-[GA]TTA[CT]GTAA[CT]-3', a sequence present in many cellular and viral promoters. Represses transcription from promoters with activating transcription factor (ATF) sites. Represses promoter activity in osteoblasts (By similarity). Represses transcriptional activity of PER1 (By similarity). Represses transcriptional activity of PER2 via the B-site on the promoter (By similarity). Activates transcription from the interleukin-3 promoter in T-cells. Competes for the same consensus-binding site with PAR DNA-binding factors (DBP, HLF and TEF) (By similarity). Component of the circadian clock that acts as a negative regulator for the circadian expression of PER2 oscillation in the cell-autonomous core clock (By similarity). Protects pro-B cells from programmed cell death (By similarity). ",
  "A - Receptor-regulated SMAD (R-SMAD) that is an intracellular signal transducer and transcriptional modulator activated by TGF-beta (transforming growth factor) and activin type 1 receptor kinases. Binds the TRE element in the promoter region of many genes that are regulated by TGF-beta and, on formation of the SMAD3/SMAD4 complex, activates transcription. Also can form a SMAD3/SMAD4/JUN/FOS complex at the AP-1/SMAD site to regulate TGF-beta-mediated transcription. Has an inhibitory effect on wound healing probably by modulating both growth and migration of primary keratinocytes and by altering the TGF-mediated chemotaxis of monocytes. This effect on wound healing appears to be hormone-sensitive. Regulator of chondrogenesis and osteogenesis and inhibits early healing of bone fractures. Positively regulates PDPK1 kinase activity by stimulating its dissociation from the 14-3-3 protein YWHAQ which acts as a negative regulator. SMAD3_HUMAN,P84022",
  "A - Binds to GC box promoters elements and selectively activates mRNA synthesis from genes that contain functional recognition sites. SP2_HUMAN,Q02086",
  "A - Transcriptional activator which forms a core component of the circadian clock. The circadian clock, an internal time-keeping system, regulates various physiological processes through the generation of approximately 24 hour circadian rhythms in gene expression, which are translated into rhythms in metabolism and behavior. It is derived from the Latin roots 'circa' (about) and 'diem' (day) and acts as an important regulator of a wide array of physiological functions including metabolism, sleep, body temperature, blood pressure, endocrine, immune, cardiovascular, and renal function. Consists of two major components: the central clock, residing in the suprachiasmatic nucleus (SCN) of the brain, and the peripheral clocks that are present in nearly every tissue and organ system. Both the central and peripheral clocks can be reset by environmental cues, also known as Zeitgebers (German for 'timegivers'). The predominant Zeitgeber for the central clock is light, which is sensed by retina and signals directly to the SCN. The central clock entrains the peripheral clocks through neuronal and hormonal signals, body temperature and feeding-related cues, aligning all clocks with the external light/dark cycle. Circadian rhythms allow an organism to achieve temporal homeostasis with its environment at the molecular level by regulating gene expression to create a peak of protein expression once every 24 hours to control when a particular physiological process is most active with respect to the solar day. Transcription and translation of core clock components (CLOCK, NPAS2, ARNTL/BMAL1, ARNTL2/BMAL2, PER1, PER2, PER3, CRY1 and CRY2) plays a critical role in rhythm generation, whereas delays imposed by post-translational modifications (PTMs) are important for determining the period (tau) of the rhythms (tau refers to the period of a rhythm and is the length, in time, of one complete cycle). A diurnal rhythm is synchronized with the day/night cycle, while the ultradian and infradian rhythms have a period shorter and longer than 24 hours, respectively. Disruptions in the circadian rhythms contribute to the pathology of cardiovascular diseases, cancer, metabolic syndromes and aging. A transcription/translation feedback loop (TTFL) forms the core of the molecular circadian clock mechanism. Transcription factors, CLOCK or NPAS2 and ARNTL/BMAL1 or ARNTL2/BMAL2, form the positive limb of the feedback loop, act in the form of a heterodimer and activate the transcription of core clock genes and clock-controlled genes (involved in key metabolic processes), harboring E-box elements (5'-CACGTG-3') within their promoters. The core clock genes: PER1/2/3 and CRY1/2 which are transcriptional repressors form the negative limb of the feedback loop and interact with the CLOCK NPAS2-ARNTL/BMAL1 ARNTL2/BMAL2 heterodimer inhibiting its activity and thereby negatively regulating their own expression. This heterodimer also activates nuclear receptors NR1D1/2 and RORA/B/G, which form a second feedback loop and which activate and repress ARNTL/BMAL1 transcription, respectively. The NPAS2-ARNTL/BMAL1 heterodimer positively regulates the expression of MAOA, F7 and LDHA and modulates the circadian rhythm of daytime contrast sensitivity by regulating the rhythmic expression of adenylate cyclase type 1 (ADCY1) in the retina. NPAS2 plays an important role in sleep homeostasis and in maintaining circadian behaviors in normal light/dark and feeding conditions and in the effective synchronization of feeding behavior with scheduled food availability. Regulates the gene transcription of key metabolic pathways in the liver and is involved in DNA damage response by regulating several cell cycle and DNA repair genes. NPAS2_HUMAN,Q99743",
  "R - Binds to a retinoid X receptor (RXR) responsive element from the cellular retinol-binding protein II promoter (CRBPII-RXRE). Inhibits the 9-cis-retinoic acid-dependent RXR alpha transcription activation of the retinoic acid responsive element. Active transcriptional corepressor of SMAD2. Links the nodal signaling pathway to the bifurcation of the forebrain and the establishment of ventral midline structures. May participate in the transmission of nuclear signals during development and in the adult, as illustrated by the down-modulation of the RXR alpha activities. TGIF1_HUMAN,Q15583",
  "A - Transcription factor that promotes adipocyte differentiation and suppresses osteoblast differentiation in the bone marrow. Enhances the osteoclast-supporting ability of stromal cells. Binds with STAT3 the consensus sequence 5'-CTTCTGGGAAGA-3' of the acute phase response element (APRE). Transactivates several promoters including FOS, OSM and PPARG. Recruits a histone deacetylase complex (By similarity).",
  "Unknown - May play an important role in B-cell differentiation as well as neural development and spermatogenesis. Involved in the regulation of the CD19 gene, a B-lymphoid-specific target gene. PAX5_HUMAN,Q02548",
  "A - Transcription factor that binds to the immunoglobulin enhancer Mu-E5/KE5-motif. Involved in the initiation of neuronal differentiation. Activates transcription by binding to the E box (5'-CANNTG-3'). Binds to the E-box present in the somatostatin receptor 2 initiator element (SSTR2-INR) to activate transcription (By similarity). Preferentially binds to either 5'-ACANNTGT-3' or 5'-CCANNTGG-3'. ITF2_HUMAN,P15884",
  "A/R - Acts as a transcriptional activator or repressor (PubMed:27181683). Plays a pivotal role in regulating lineage-specific hematopoiesis by repressing ETS1-mediated transcription of erythroid-specific genes in myeloid cells. Required for monocytic, macrophage, osteoclast, podocyte and islet beta cell differentiation. Involved in renal tubule survival and F4/80 maturation. Activates the insulin and glucagon promoters. Together with PAX6, transactivates weakly the glucagon gene promoter through the G1 element. SUMO modification controls its transcriptional activity and ability to specify macrophage fate. Binds element G1 on the glucagon promoter (By similarity). Involved either as an oncogene or as a tumor suppressor, depending on the cell context. MAFB_HUMAN,Q9Y5Q3",
  "A - Transcription activator that binds DNA cooperatively with DP proteins through the E2 recognition site, 5'-TTTC[CG]CGC-3' found in the promoter region of a number of genes whose products are involved in cell cycle regulation or in DNA replication. The DRTF1/E2F complex functions in the control of cell-cycle progression from G1 to S phase. E2F4 binds with high affinity to RBL1 and RBL2. In some instances can also bind RB1. Specifically required for multiciliate cell differentiation: together with MCIDAS and E2F5, binds and activate genes required for centriole biogenesis.")
functions <- transpose(as.data.frame(functions_u))

cis_gene_cards_fxns <- data.frame(z, 
                             unique(cis_gene_cards),
                             functions)
colnames(cis_gene_cards_fxns)[6:7] <- c("Link", "Function")
cis_gene_cards_fxns$Function_Abbreviation <- gsub(" \\-.*$", "", cis_gene_cards_fxns$Function)
cis_gene_cards_fxns$CisZSign <- sign((cis_gene_cards_fxns$cisZscore))
cis_gene_cards_fxns$TransZSign <- sign((cis_gene_cards_fxns$transZscore))

# Expected:
# TF    cis-G   trans-G
# Act   z > 0   z > 0
# Rep   z < 0   z > 0
# Rep   z > 0   z < 0
# Act   z < 0   z < 0

# Does not take into account A/R
x <- cis_gene_cards_fxns

#backup: 
cis_gene_cards_fxns <- x
cis_gene_cards_fxns <- merge(x, data.frame(snps_interesting_TFmatch[,"SNP"], snps_interesting_TFmatch[,"cisGene"]),by.x = "cisGene", by.y="cisGene", all.y=T)
cis_gene_cards_fxns <- as.data.frame(cis_gene_cards_fxns)
for (i in 1:length(cis_gene_cards_fxns)){
  if (cis_gene_cards_fxns$CisZSign[i]*cis_gene_cards_fxns$TransZSign[i] == 1)
    cis_gene_cards_fxns$Calculated_Function[i]<- "A" 
  else 
    cis_gene_cards_fxns$Calculated_Function[i]  <- "R"
}
# cis_gene_cards_fxns$Calculated_Function <- apply(cis_gene_cards_fxns, 1, FUN=function(x) if (x$CisZSign*x$TransZSign == 1) "A" else "R")
for (i in 1:dim(cis_gene_cards_fxns)[1]){
  if (identical(cis_gene_cards_fxns$Function_Abbreviation[i],cis_gene_cards_fxns$Calculated_Function[i]))
  {
    cis_gene_cards_fxns$DoesItMatch[i] <- 1 #T
    #print(paste(cis_gene_cards_fxns$Function_Abbreviation[i], " ", cis_gene_cards_fxns$Calculated_Function[i], " ", cis_gene_cards_fxns$DoesItMatch[i]), " paste 1")
  }
  else 
  {
    cis_gene_cards_fxns$DoesItMatch[i]  <- 0 #F
    #print(paste(cis_gene_cards_fxns$Function_Abbreviation[i], " ", cis_gene_cards_fxns$Calculated_Function[i], " ", cis_gene_cards_fxns$DoesItMatch[i]), " paste 0")
  }
}

# NOT COMPLETELY ACCURATE BECAUSE A/R IS UNKNOWN
cis_gene_cards_fxns %>% count(DoesItMatch)
# 0 = F = Doesn't Match  = 358
# 1 = T = Does Match = 100

y <- aggregate(cis_gene_cards_fxns$DoesItMatch, by=list(SNP=cis_gene_cards_fxns$SNP), FUN=mean)

# Accuracy of Prediction
# SNP         x
# 1    rs10445308 1.0000000
# 2     rs1055348 0.0000000
# 3    rs10744777 1.0000000
# 4    rs10774624 0.3750000
# 5    rs10774625 0.4000000
# 6    rs10821936 1.0000000
# 7    rs10844706 0.0000000
# 8    rs11052877 0.0000000
# 9    rs11065979 0.4000000
# 10   rs11065987 0.4000000
# 11   rs11066301 0.2500000
# 12   rs11078927 1.0000000
# 13  rs111230933 0.6666667
# 14  rs111255518 0.0000000
# 15  rs111410428 0.0000000
# 16  rs111496944 0.0000000
# 17   rs11204677 0.0000000
# 18  rs112115374 0.0000000
# 19  rs112227868 0.0000000
# 20  rs112296382 0.0000000
# 21  rs112954547 0.0000000
# 22  rs113589048 0.0000000
# 23  rs114012716 0.0000000
# 24  rs114127008 0.0000000
# 25  rs114211054 0.0000000
# 26  rs114240154 0.0000000
# 27  rs114250120 0.5000000
# 28  rs114260710 0.0000000
# 29  rs114348871 0.0000000
# 30  rs114394153 0.0000000
# 31  rs114415823 0.0000000
# 32  rs114487324 0.0000000
# 33  rs114535641 0.0000000
# 34  rs114571307 0.0000000
# 35  rs114573742 0.0000000
# 36  rs114592706 0.0000000
# 37  rs114638960 0.0000000
# 38  rs114874012 0.0000000
# 39  rs114962780 0.0000000
# 40  rs114987030 0.0000000
# 41  rs115029884 0.0000000
# 42  rs115054796 0.0000000
# 43    rs1150753 0.0000000
# 44    rs1150757 0.0000000
# 45  rs115116967 0.0000000
# 46  rs115122307 0.2500000
# 47  rs115146037 0.0000000
# 48  rs115164593 0.0000000
# 49  rs115318382 0.0000000
# 50  rs115326065 0.0000000
# 51  rs115344853 0.0000000
# 52  rs115347165 0.0000000
# 53  rs115374828 0.0000000
# 54  rs115484360 0.0000000
# 55  rs115509556 0.3333333
# 56  rs115603641 0.0000000
# 57  rs115610190 0.0000000
# 58  rs115611791 0.0000000
# 59  rs115613956 0.0000000
# 60  rs115619714 0.0000000
# 61  rs115648972 0.0000000
# 62  rs115687010 0.2500000
# 63  rs115715453 0.0000000
# 64  rs115741842 0.0000000
# 65  rs115772457 0.0000000
# 66  rs115845232 0.2500000
# 67  rs115902351 0.0000000
# 68  rs116026314 0.0000000
# 69  rs116152465 0.0000000
# 70  rs116212130 0.0000000
# 71  rs116336456 0.0000000
# 72  rs116351884 0.0000000
# 73  rs116358946 0.0000000
# 74  rs116392568 0.0000000
# 75  rs116396237 0.0000000
# 76  rs116414615 0.0000000
# 77  rs116528482 0.0000000
# 78  rs116627865 0.0000000
# 79  rs116628560 0.0000000
# 80  rs116633882 0.0000000
# 81  rs116743258 0.0000000
# 82  rs116766239 0.0000000
# 83  rs116766442 0.0000000
# 84  rs116778584 0.0000000
# 85  rs116839689 0.0000000
# 86   rs11696739 0.0000000
# 87   rs12485738 0.0000000
# 88    rs1265097 0.7500000
# 89    rs1265564 0.5000000
# 90    rs1270942 0.0000000
# 91   rs12946510 1.0000000
# 92   rs13015714 1.0000000
# 93   rs13325613 0.0000000
# 94    rs1354034 0.0000000
# 95    rs1420103 1.0000000
# 96  rs144033971 0.0000000
# 97  rs149110519 0.0000000
# 98  rs150881176 1.0000000
# 99    rs1569723 0.0000000
# 100    rs174555 1.0000000
# 101  rs17613465 0.0000000
# 102  rs17696736 0.3750000
# 103   rs1799964 0.0000000
# 104 rs182050989 0.0000000
# 105   rs1883832 0.0000000
# 106   rs1892548 0.0000000
# 107   rs2058660 1.0000000
# 108   rs2071591 0.0000000
# 109   rs2229094 0.0000000
# 110   rs2305480 1.0000000
# 111    rs234709 1.0000000
# 112   rs2523607 0.0000000
# 113   rs2596500 0.0000000
# 114   rs2617170 0.1666667
# 115  rs28453840 0.0000000
# 116   rs2872507 1.0000000
# 117  rs28798705 0.0000000
# 118   rs3130564 1.0000000
# 119   rs3130614 0.0000000
# 120   rs3131379 0.0000000
# 121   rs3184504 0.3846154
# 122   rs3774937 0.0000000
# 123   rs3809272 0.5000000
# 124   rs3809627 0.0000000
# 125   rs3819306 0.0000000
# 126   rs4239702 0.0000000
# 127   rs4328821 1.0000000
# 128   rs4766578 0.4166667
# 129    rs477515 0.0000000
# 130   rs4785587 0.0000000
# 131   rs4794820 0.5000000
# 132   rs4795397 1.0000000
# 133   rs4810485 0.0000000
# 134    rs492400 0.0000000
# 135   rs4947311 0.0000000
# 136   rs4970966 0.0000000
# 137  rs59716545 1.0000000
# 138    rs597808 0.3636364
# 139   rs6032662 0.0000000
# 140   rs6074022 0.0000000
# 141    rs615672 0.0000000
# 142   rs6419573 1.0000000
# 143    rs653178 0.3636364
# 144   rs6563831 0.0000000
# 145   rs6586282 1.0000000
# 146   rs6695223 0.0000000
# 147   rs6708413 1.0000000
# 148    rs674313 0.0000000
# 149   rs6785206 1.0000000
# 150   rs6929796 0.0000000
# 151   rs7089424 1.0000000
# 152   rs7090445 1.0000000
# 153   rs7210990 0.0000000
# 154   rs7310615 0.3846154
# 155  rs74942078 0.0000000
# 156   rs7570971 0.0000000
# 157   rs7616215 0.0000000
# 158  rs77216612 1.0000000
# 159   rs7740107 0.0000000
# 160   rs7776054 1.0000000
# 161   rs7973618 0.0000000
# 162   rs8067378 1.0000000
# 163   rs8069176 1.0000000
# 164    rs917997 1.0000000
# 165   rs9273076 0.0000000
# 166   rs9376090 1.0000000
# 167   rs9399137 1.0000000
# 168    rs990171 1.0000000

# y %>% count(x)
# x     n
# <dbl> <int>
# 1 0       115
# 2 0.167     1
# 3 0.25      4
# 4 0.333     1
# 5 0.364     2
# 6 0.375     2
# 7 0.385     2
# 8 0.4       3
# 9 0.417     1
# 10 0.5       4
# 11 0.667     1
# 12 0.75      1
# 13 1        31

# y <- aggregate(cis_gene_cards_fxns$DoesItMatch, by=list(Gene=cis_gene_cards_fxns$cisGene_commonname), FUN=mean)
# Gene x
# 1   BCL11A 0
# 2  BHLHE40 0
# 3     E2F4 0
# 4     EBF1 0
# 5     EGR2 0
# 6    EOMES 0
# 7    FOXO1 0
# 8     IRF4 0
# 9     KLF3 0
# 10    KLF4 0
# 11    KLF6 1
# 12    MAFB 0
# 13     MAX 0
# 14   MEF2C 1
# 15   MEIS1 0
# 16     MYB 1
# 17   NFIL3 1
# 18   NPAS2 1
# 19    PAX5 0
# 20    RARA 0
# 21    RORA 0
# 22   SMAD3 0
# 23     SP2 1
# 24   STAT1 1
# 25    TCF4 0
# 26   TGIF1 0
# 27    ZEB1 0
# 28  ZNF467 0

# y %>% count(x)
# A tibble: 2 x 2
# x     n
# <dbl> <int>
# 1     0    21
# 2     1     7

# ==========================================================
# Trans GeneCards
# ==========================================================

trans_gene_cards <- data.table(snps_interesting_TFmatch[,unique(snps_interesting_TFmatch$transGene)], paste("https://www.genecards.org/cgi-bin/carddisp.pl?gene=", trim(snps_interesting_TFmatch[,unique(snps_interesting_TFmatch$transGene)]), sep=""))

# Install biomaRt
# source("https://bioconductor.org/biocLite.R")
# biocLite("biomaRt")

# https://www.biostars.org/p/178726/
library("biomaRt")

mart <- useMart("ENSEMBL_MART_ENSEMBL")
mart <- useDataset("hsapiens_gene_ensembl", mart)
ens <- snps_interesting_TFmatch[,unique(snps_interesting_TFmatch$transGene)]
ensLookup <- gsub("\\.[0-9]*$", "", ens) # not necessary

annotLookup <- getBM(
  mart=mart,
  attributes=c("ensembl_transcript_id", "ensembl_gene_id", "gene_biotype", "external_gene_name"),
  filter="ensembl_gene_id",
  values=ensLookup,
  uniqueRows=TRUE)

annotLookup2 <- data.frame(
  ens[match(annotLookup$ensembl_gene_id, ensLookup)],
  annotLookup)

colnames(annotLookup2) <- c(
  "original_id",
  c("ensembl_transcript_id", "ensembl_gene_id", "gene_biotype", "external_gene_name"))
annotLookup_abbrev <- data.table(unique(annotLookup2[, "original_id"]), unique(annotLookup2[, "external_gene_name"]))
fwrite(annotLookup_abbrev, "q2_1_written/trans_gene_commonname.csv")