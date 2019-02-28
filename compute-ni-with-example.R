# Copyright 2017-2018. Triad National Security, LLC. All rights reserved.

require(glm2)
require(beeswarm)
require(RColorBrewer)

compute.ni <- function(X, Y=NULL, serum.name=NULL, virus.names=NULL, invert=T) {

    # The invert flag is used for handling both sera and antibodies.
    # Use invert=T for sera, where increasing titer means decreased potency.
    # Use invert=F for antibodies, where increasing titer means increased potency.
    # It is kind of a pain but enables using the same criteria and plotting code.

    if (is.matrix(X)) {
        if (ncol(X) == 2) {
            if (is.null(serum.name))
                serum.name = colnames(X)[1]
            if (is.null(virus.names))
                virus.names = rownames(X)
            if (!is.null(Y))
                stop("Over-writing Y value because X has two columns.")
            Y = X[, 2]
            X = X[, 1]
        } else if (ncol(X) == 1 & !is.null(Y)) {
            if (is.null(serum.name))
                serum.name = colnames(X)[1]
            if (is.null(virus.names))
                virus.names = rownames(X)
        } else if (nrow(X) == 1 & !is.null(Y)) {
            if (is.null(serum.name))
                serum.name = rownames(X)[1]
            if (is.null(virus.names))
                virus.names = colnames(X)
            X = t(X)
        } else { # X is a matrix with more than 2 columns or 1 row
            stop("ERROR: not sure how to interpret the matrix passed to compute.ni()")
        }
    }
    if (length(X) != length(Y))
        stop("ERROR: The arguments passed to compute.ni() are of unequal lengths")
    if (all(is.na(X)) | all(is.na(Y)) | all(is.null(X)) | all(is.null(Y)))
        stop("ERROR: Arguments passed to compute.ni() contain no data")
    if (!is.logical(X) | !is.numeric(Y))
        stop("ERROR: Arguments passed to compute.ni() are incompatible with expectations")
    # verify that X is logical
    # verify that Y is numeric
    ### invert == F, is_neutralized -> T, F ~ 0, 5  id50 > cutoff
    ### invert == T, is_resistant   -> F, T ~ 0, 5  id50 < cutoff

    if (invert)
        X = !X

    this.fit <- try(glm(c(!invert, invert, X) ~ c(0, 5, Y), family='binomial'))
    # this.fit <- try(glm(c(T, F, X) ~ c(0, 5, Y), family='binomial'))

    if (!this.fit$converged)
        this.fit <- try(glm2::glm2(c(!invert, invert, X) ~ c(0, 5, Y),
            binomial, control=list(maxit=500)))

    this.ni = list()

    test.c = try(anova(this.fit, test='Chisq'))
   # test.c = try(anova(this.fit, test='Rao'))


    if (is.null(virus.names) & !is.null(names(Y)))
        virus.names = names(Y)

    if (length(test.c$"Pr(>Chi)") == 2)
        this.ni$P = test.c$"Pr(>Chi)"[[2]]#/2

    this.ni$is.inverted = invert
    this.ni$X = X
    this.ni$Y = Y
    this.ni$serum.name = serum.name
    this.ni$virus.names = virus.names

    this.ni$index = -this.fit$coefficients[[1]] / this.fit$coefficients[[2]]
    this.ni$beta0 = this.fit$coefficients[[1]]
    this.ni$beta1 = this.fit$coefficients[[2]]
    this.ni$deviance = this.fit$deviance
    this.ni$null.deviance = this.fit$null.deviance

    class(this.ni) = "ni"

    return(this.ni)
}

# make beeswarm panels for the selected serum
plot.ni <- function(this.ni, show.virus.names=T, ...) {

    if (class(this.ni) != "ni")
        return (NULL)
    
    dots = list(...)

    if (!is.null(dots$xlab)) { x.lab = dots$xlab } else {x.lab = 'Tiered Virus Geometric Mean ID50'}
    if (!is.null(dots$xlim)) { my.xlim = dots$xlim } else { my.xlim <- c(1/2,3.75) }

    ### NB: tiered.virus.means is defined outside of this environment
    virus_colors <- colorRampPalette(c("grey",
        RColorBrewer::brewer.pal(9, "Reds")[3:9]))(n = length(tiered.virus.means))
    names(virus_colors) = names(tiered.virus.means)

    serum_colors <- colorRampPalette(c("grey",
                    RColorBrewer::brewer.pal(9, "Blues")[3:9]))(n = 101)
    names(serum_colors) = sprintf("%4.2f", 5*c(0:100)/100)
    virus.colors = "black"
    if (all(this.ni$virus.names %in% names(virus_colors))) {
        virus.colors = virus_colors[this.ni$virus.names]
    } else {
        virus.colors = virus_colors[colnames(this.ni$X)]
    }

    ### ID50 is correct for sera but not for bnAbs:
    y.lab = ifelse(this.ni$is.inverted,
                   paste0('ID50 <', my.cutoff, '?'),  # is_resistant
                   paste0('ID50 >', my.cutoff, '?') ) # is_neutralized
#    y.lab = ifelse(this.ni$is.inverted,
#                   "Is Env Resistant?",  # is_resistant
#                   "Is Env Neutralized?") #) # is_neutralized

        X.ordered = ordered(this.ni$X, c(F, T))
#    if (this.ni$is.inverted)
#        X.ordered = ordered(this.ni$X, c(T, F))

    bees.out <- beeswarm::beeswarm(
        this.ni$Y ~ X.ordered,
        pwcol=virus.colors,
        pch=20, #cex=2/3,
        xlim=my.xlim,
        #labels=NA,
        horizontal=T, method='hex',	corral='wrap',
        ylab=y.lab,
        xlab=x.lab, bty='n')

    serum_color = serum_colors[as.character(sprintf("%4.2f", round((this.ni$index*100)/5)/20))]
    abline(v=this.ni$index, col=serum_color, lty=2, lwd=3/2)

    n.points <- 100

    lines(c((my.xlim[1]*n.points):(my.xlim[2]*n.points))/n.points,
            2 - 1/(1+exp( (this.ni$beta0 +
            # in the beeswarm view, y=0 becomes y=1, so shift curve up that amount
            this.ni$beta1 * c((my.xlim[1]*n.points):(my.xlim[2]*n.points))/n.points))),
              col=serum_color, lwd=3/2)

        abline(h=3/2, col='#77777777')

    usr=par('usr')
    y = c(1,2)
    if (!this.ni$is.inverted)
        y = rev(y)

    my.cex.lab=par('cex.lab')
    mtext("Resistant", 2, at=y[2], line=-1/2, font=2, cex=my.cex.lab, padj=0)
    mtext("Neutralized", 2, at=y[1], line=-1/2, font=2, cex=my.cex.lab, padj=0)

    if (show.virus.names & !all(is.null(this.ni$virus.names)) &
        length(this.ni$virus.names) == nrow(bees.out)) {

        for (i in 1:nrow(bees.out)) {

            this.adj = 0

            # attempt to reduce overplotting among virus names
            # we consider the length of the x.orig vector to handle cases where all assays
            # are positive and thus x values are centered on 1 not 2
            if ((bees.out$x.orig[i] == T & bees.out$x[i] < length(unique(bees.out$x.orig))) |
                (bees.out$x.orig[i] == F & bees.out$x[i] < 1))
                this.adj = 1

            text(bees.out$y[i], bees.out$x[i],
                 paste0("  ", this.ni$virus.names[which(this.ni$Y == bees.out$y.orig[i])], "  "),
                 srt=60, adj=this.adj, cex=1/2)
        }
    }

    if (!is.null(this.ni$serum.name))
        mtext(paste0('  ',
                this.ni$serum.name),
                side=3, adj=0, line=-5/4)#, cex=4/5)
    if (!is.na(this.ni$index))
        mtext(paste0('NI=',
            signif(this.ni$index, 3), " "), side=3, adj=1, line=-5/4)#, cex=4/5)
#    if (!is.na(this.ni$beta0))
#        mtext(paste0(' Intercept=',
#            signif(this.ni$beta0, 3)), side=1, adj=0, line=-5/4)#, cex=4/5)
#    if (!is.na(this.ni$beta1))
#        mtext(paste0('Slope=',
#            signif(this.ni$beta1, 3), " "), side=1, adj=1, line=-5/4)#, cex=4/5)
    usr=par('usr')
    xpos = ifelse(this.ni$index <= mean(usr[1:2]),
                  ifelse(this.ni$index < usr[1], usr[1], this.ni$index),
                  ifelse(this.ni$index > usr[2], usr[2], this.ni$index))
    ypos = ifelse(this.ni$index <= mean(usr[1:2]), 1.005*3/2, 2-1.005*3/2)
    xadj1 = ifelse(this.ni$index <= mean(usr[1:2]), 0, 1)
    xadj2 = ifelse(this.ni$index <= mean(usr[1:2]), 1, 0)

    if (!is.null(this.ni$P))
        text(xpos, 1.005*3/2, paste0(' P=', signif(this.ni$P, 3)), 
             adj=c(xadj1, xadj2))

    return (bees.out)
}


if (!exists("tiered.virus.means")) {

    tiered.virus.means = c(3.4518952, 3.4500124, 3.4018602,
    3.3931649, 3.3917412, 3.3818863, 3.3777216, 3.3431046,
    3.3274112, 3.3249167, 3.2981185, 3.2894911, 3.2877544,
    3.2870648, 3.274706, 3.2719438, 3.2695783, 3.2382615,
    3.2349504, 3.2270931, 3.210122, 3.2068208, 3.1942505, 
    3.18653, 3.1864724, 3.1855454, 3.1752743, 3.169986, 
    3.1618635, 3.1518948, 3.1513432, 3.1433353, 3.1425163, 
    3.1413561, 3.1358956, 3.1211911, 3.1118822, 3.0977091, 
    3.0961647, 3.0944829, 3.0893667, 3.0815745, 3.0810291, 
    3.0746799, 3.0732143, 3.0651357, 3.0633387, 3.060129, 
    3.0578574, 3.0565372, 3.0551376, 3.0527407, 3.0495026, 
    3.0484212, 3.0365098, 3.0343924, 3.0327621, 3.031708, 
    3.0287945, 3.025137, 3.0246854, 3.02095, 3.0169311, 
    3.0140457, 3.0122815, 3.0115831, 3.0046311, 2.9982194, 
    2.9832534, 2.983044, 2.9809387, 2.9737042, 2.9736247, 
    2.9672707, 2.9631549, 2.9622094, 2.9589299, 2.9573763, 
    2.9557248, 2.9551372, 2.9522923, 2.9504675, 2.9474156, 
    2.9455703, 2.9432987, 2.9310755, 2.9307542, 2.9284452, 
    2.9273051, 2.924236, 2.9236456, 2.923552, 2.9153241, 
    2.9120738, 2.9037101, 2.9018558, 2.9000251, 2.8997145, 
    2.8969792, 2.8943559, 2.8932564, 2.8907908, 2.8874433, 
    2.8864076, 2.8832747, 2.881502, 2.8785686, 2.8781595, 
    2.8666252, 2.8591159, 2.8527228, 2.8443209, 2.8407145, 
    2.8400839, 2.8311313, 2.8300269, 2.8176766, 2.8130717, 
    2.8128497, 2.8128387, 2.8064153, 2.8057773, 2.8015693, 
    2.7995136, 2.7962432, 2.7850592, 2.7825023, 2.7821747, 
    2.7809088, 2.7786472, 2.7734253, 2.7725501, 2.769309, 
    2.7648532, 2.7624087, 2.7593683, 2.7509074, 2.7454313, 
    2.7434308, 2.7431746, 2.7170641, 2.7090047, 2.7078157, 
    2.6949421, 2.6946144, 2.6874032, 2.6855495, 2.6772897, 
    2.6654513, 2.649899, 2.6458381, 2.6441812, 2.641831, 
    2.6411945, 2.6396375, 2.6394637, 2.63484, 2.6336714, 
    2.6329086, 2.6255975, 2.6187343, 2.6123852, 2.6123832, 
    2.6104459, 2.6096795, 2.6088344, 2.6056025, 2.5993445, 
    2.5972838, 2.595923, 2.5953628, 2.5939003, 2.5937991, 
    2.593061, 2.5886124, 2.5883879, 2.5847263, 2.5758654, 
    2.5687021, 2.5674827, 2.5661493, 2.5629708, 2.5595255,
    2.5562236, 2.5545878, 2.5498402, 2.5490415, 2.527825, 
    2.5165746, 2.5113306, 2.5012773, 2.4991492, 2.4950373,
    2.4833008, 2.4805315, 2.4668383, 2.4474403, 2.4299406,
    2.4284314, 2.4224879, 2.4147885, 2.4007442, 2.3966842,
    2.3804308, 2.3800018, 2.3658578, 2.364776, 2.3217483,
    2.2970873, 2.2927141, 2.2766859, 2.2696272, 2.2693301,
    2.212634, 2.1938917, 2.1596216, 2.0868452, 2.0446179, 
    2.0415727, 2.0401989, 1.9404463, 1.8006818, 1.6620417,
    1.6574762, 1.1988803)

    names(tiered.virus.means) = c( "9014_01_TB1_4769", "98-F4_H5-13",
 	"231965.C01", "QH0515.1", "569-F1_37_10", "6471.V1.C16",
 	"3637.V5.C3", "X2088_C9", "PSR-0508.2", "X2160_C25",
 	"T251-18", "CE2010_F5", "CE2103_E8", "NKR-0512.8", "CNE28",
 	"3016.V5.C45", "63358_04_P3_4013", "BG1168.1", "541-F1_A7_2",
 	"CE706010018_2E3", "3168.V4.C10", "1051_12_C22_3325",
 	"192018_B1_9", "6240_08_TA5_4622", "X1254_C3", "Q461.E2",
 	"CNE30", "193003_B10", "1054_07_TC4_1499", "6022.V7.C24",
 	"62357_14_D3_4589", "CNE5", "KSS-0514.13", "A07412M1.VRC12",
 	"Q769.D22", "CH120.6", "6101.10", "BJOX025000.01.1", "CNE23",
 	"RPW-0510.2", "7165.18", "TRJO4551.58", "X1100_C7",
 	"0815.V3.C3", "410-F2_1_30", "3817.V2.C59", "BF1677F2.613A",
 	"CE703010054_2A2", "BORI_D9_4D7_1410", "ZM214M.PL15",
 	"H086.8", "H035.18", "CNE15", "700010040_C9_4520", "PVO.4",
 	"CNE55", "1012_11_TC21_3257", "ZM53M.PB12", "620345.C01",
 	"T257-31", "TT31P_2F10_2792", "191084_B7-19", "Q259.D2.17",
 	"CNE6", "1058_11_B11_1550", "HIV-16055-2.3", "H078.14",
 	"191727_F6", "R2184.C04", "H080.23", "CNE52",
 	"CE703010217_B6", "3326.V4.C3", "1051_12_TD12_3291",
 	"WITO4160.33", "CNE31", "ZM247V1(REV-)", "1656.P21", "H030.7",
 	"252-7", "21020_13", "BORI_D9_4F8_1413", "1006_11_C3_1601",
 	"3728.V2.C6", "6545.V4.C1", "CE704010083_B8", "235-47",
 	"269-12", "191821_E6_1", "C2101.C01", "HIV-26191-2.48",
 	"249M_B10", "CE0682_E4", "CH115.12", "T278-50", "SC422661.8",
 	"TT29P_3A1_2769", "CAP210.2.00.E8", "BF1266.431A",
 	"AC10.0.29", "HIV-25711-2.4", "6040.V4.C15",
 	"BJOX009000.02.4", "SC45_4B5_2631", "CH070.1", "T266-60",
 	"3468.V1.C12", "CNE18", "CNE11", "3103.V3.C10",
 	"CAP45.2.00.G3", "X2252_C7", "6244_13_B5_4576", "X1193_C1",
 	"CE2060_G9", "C4118.C09", "6811.V7.C18", "816763.C02",
 	"3301.V1.C24", "PRB931_06_TC3_4930", "SPK-0525.13",
 	"Q842.D12", "HIV-16936-2.21", "CNE21", "1105.P17.1",
 	"6480.V4.C25", "HIV-00836-2.5", "246F_C1G", "QH0692.42",
 	"CAAN5342.A2", "CNE17", "X1854_C2_10", "6785.V5.C14",
 	"477-F3_13_55", "BJOX019000.02.1", "BJOX028000.10.3",
 	"402288.C01", "BJOX015000.11.5", "CE1172_H1",
 	"1059_09_A4_1460", "CE0393_C3", "DU422.1",
 	"700010058_A4_4375", "6041.V3.C23", "1056_10_TA11_1826",
 	"T280-5", "ZM249M.PL1", "HIV-16845-2.22", "R1166.C01",
 	"THRO4156.18", "7060101641A7(REV-)", "BF942.218D",
 	"191955_A11", "CE703010228_1C4", "CNE16", "HIV-25925-2.22",
 	"263-8", "6952.V1.C20", "BJOX010000.06.2", "ZM233M.PB6",
 	"CH038.12", "X1632_S2_B10", "211-9", "427299.C12", "CH119.10",
 	"CE704010042_2E5", "6540.V4.C1", "9004SS_A3_4", "3169.P4",
 	"H061.14", "MS208.A1", "DU172.17", "CH181.12", "T255-34",
 	"ZM135M.PL10A", "270015__J5_1", "HIV-001428-2.42",
 	"RHPA4259.7", "C1080.C03", "7030102001E5(REV-)", "2705.P18.1",
 	"3273.P21.1", "246-F3_C10_2", "CE1176_A3",
 	"PRB958_06_TB1_4305", "SC05_8C11_2344", "DU151.2", "21023_B4",
 	"PWJ-0513.39", "R18553_E1", "CH117.4", "CNE9", "234-F1-16-57",
 	"DU123.6", "9020_20_A13_4607", "TRO.11", "BJOX002000.03.2",
 	"H029.12", "P1981_C5_3", "X2278_C2_B6", "C3347.C11",
 	"WEAU_D15_410_5017", "CH111.8", "CE704809221_1B3",
 	"CE1086_B2", "401-F1_8_10", "REJO4541.67", "TV1.21",
 	"CH114.8", "6535.3", "T250-4", "242-14", "ZM109F.PB4", "CNE8",
 	"ZM197M.PB7", "1394C9G1(REV-)", "CNE19", "HIV-25710-2.43",
 	"H031.7", "DU156.12", "CNE20", "398-F1_F6_20", "92BR025.9",
 	"CNE67", "6644.V2.C33"
 	)
}

show.example = T

if (show.example) {
    example.in <- read.csv("c_clade_sample.csv", header=T, stringsAsFactors=F)

    # clean up column names corrupted by R - 
    colnames(example.in) = gsub("^X", "", gsub(".REV..", "(REV.)", colnames(example.in), fixed=T))
    # match names of tier-scored Envs to input file
    tmp.env.names = gsub('-', '.', names(tiered.virus.means), fixed=T)
    matched.envs = sapply(1:ncol(example.in), function(i)
        ifelse(colnames(example.in)[i] %in% tmp.env.names, 
               which(tmp.env.names == colnames(example.in)[i]), NA))

    if (any(is.na(matched.envs)))
        stop(paste("Unable to find an Env to match", colnames(example.in)[which(is.na(matched.envs))], "\n"))

    Y <- tiered.virus.means[matched.envs] # these (viruses) are the same for all rows

    my.cutoff=50 # ID50s below this value are considered "not neutralized"

    # add column/s to the input data matrix
    example.in["nis"] <- NA
    example.in$nis <- sapply(1:nrow(example.in), function(i) {
        X <- unlist(example.in[i, which(!is.na(matched.envs))] > my.cutoff)
        serum_name = rownames(example.in)[i]
        (retval <- compute.ni(X, Y, serum.name=serum_name))$index 
        } )

    # randomly choose one serum and plot it, as an example
    plot.example = sample(nrow(example.in), 1)
    X = as.numeric(example.in[plot.example, -ncol(example.in)])
    names(X) = names(Y)
    this.result = compute.ni(cbind(t(X > my.cutoff)), Y)
    this.plot = plot.ni(this.result)
}
