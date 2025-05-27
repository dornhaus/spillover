# Alasdair Houston & Anna Dornhaus
# R script for numerical model analysis and figures
# Spillover project: Bees 'spilling over' into crop from wildflowers planted to attract them

# Paper reference: XXX (to be updated on acceptance of the paper)
# PAPER VERSION R1 April 2025
# SUPPLEMENTARY MATERIALS

### Contents ------------------
# I Models A.1, A.2, A.3
# II Model B

### Libraries and graphics setup ------------------
# Libraries purely for graphics
library(scales)
library(viridis)
### COLORS 1
# colorlist has the colors for 
# N, Fc, Nc, D, Fw, Nw, the tangent in this order
templist1 <- inferno(6)
templist2 <- mako(6)
colorlist <- c(templist1[2], templist1[3], templist1[4], templist1[5], templist2[4], templist2[5], templist2[3])

### COLORS 2
# colors for parameter variations
no_states <- 5
par_colors_dark <- viridis(no_states)
par_colors <- alpha(par_colors_dark, 0.5)
N_colors <- mako(no_states, begin = 0.3, end = 0.9)
F_colors <- magma(no_states, begin = 0.3, end = 0.9)

#dev.off()
opar <- par()


### PARAMETERS ------------------

# Default parameter set
c_fixed <- 8

aw <- 0.4 #diminishing returns of attraction based on 'amount' of flowers 
AAw <- 10   # constant translating 'amount' of wildflowers into pollinators
AAc <- 1 # constant translating 'amount' of crop into pollinators
ac <- 0.4
l <- 10
rho <- 0.02
gamma <- 3
 
## For numerical calculations
max_w <- 1000 # for numerical calculations; min >2
resolution <- 10000 # number of calculated values; min 200
w_list <- c(seq(from=0, to=2, length.out=200), seq(from=2, to=max_w, length.out=resolution-200))
c_list <- w_list

## Parameter variations
# no_states needs to be defined before colors
c_states <- round(seq(1, 20, length.out=no_states), 1)
w_states <- c_states
aw_list <- round(seq(0.01,1, length.out=no_states), 2) 
ac_list <- round(seq(0.01,1, length.out=no_states), 2) 
range_param <- 1
expon <- seq(-range_param, range_param, length.out=no_states)
AAw_list <- round(10^expon, 2)
AAc_list <- round(10^expon, 2) 

### MODELS A ------------------
### GENERAL FUNCTIONS -----------------
# Equation for D
D <- function(c,w) {
  return(c/(c+w))
}
# Equation for F
F <- function(N, D) {
  return(N*D)
}
##### MODEL A.1 ----------------------------------
N_A1 <- function(c,w) {
  return((AAw*w) / (aw+w) + (AAc*c) / (ac+c))
}

##### MODEL A.2 ----------------------------------
N_A2 <- function(c,w) {
  return(AAw*(w^aw) + AAc*(c^ac))
}

##### MODEL A.3 ----------------------------------
N_A3 <- function(c,w) {
  return(AAw*(1-exp(-aw*w)) + AAc*(1-exp(-ac*c)))
}



### NUMERICAL CALCULATIONS ------------------------

# First, we calculate outcomes across w for a fixed c
# (nothing optimized here, no w*)
D_over_w <- D(c_fixed, w_list)
Ns_A1 <- N_A1(c_fixed, w_list)
Ns_A2 <- N_A2(c_fixed, w_list)
Ns_A3 <- N_A3(c_fixed, w_list)
Fs_A1 <- F(Ns_A1, D_over_w)
Fs_A2 <- F(Ns_A2, D_over_w)
Fs_A3 <- F(Ns_A3, D_over_w)

### FIG S1 ATTRACTION VS CHOICE OVER D - 3X3 ------------------------
# Figure format -------------------------------
## For plotting
pointsize <- 0.6
Nlim <- 25 # Max for N graph
F_clim <- Nlim / 10 # Max for F/c graph
wlim <- 20 # for plotting
clim <- 50 # Max c for graph
offsetx <- 0.1
offsety <- 0.9
point_transparency <- 0.5

# Make 3x3 panels
par(mfrow=c(3,3))
# Various other format adjustments (margins in order bottom, left, top, right)
par(oma = c(5,5,1,0), mar=c(0.5, 0.5, 0.5, 0.5), mgp=c(3, 1, 0), las=1)

# 1st ROW - Number of bees, F, Nc, w* -----------------

# ROW 1: Number of bees: N, Fc, Fw=N-Fc, Nc
# column 1: Model A1
plot(NULL
     , ylim= c(0, Nlim)
     , xlim= c(0, 1)
     , xaxt= 'n'
)
# Fc filled area
polygon(c(D_over_w, rev(D_over_w)), c(c(1:resolution)*0, rev(Fs_A1))
        , col = alpha(colorlist[2], point_transparency)
        , lty = 0
)
# Fw = area between Fc and N
polygon(c(D_over_w, rev(D_over_w)), c(Fs_A1, rev(Ns_A1))
        , col = colorlist[5]
        , lty = 0
)
# horizontal line at Nc
abline(h=N_A1(c_fixed, 0), lty=2, lwd=3, col=colorlist[3])
# Cross of lines at w*
abline(h=max(Fs_A1), lwd=2, lty=3, col="grey")
abline(v=D_over_w[which(Fs_A1==max(Fs_A1))], lwd=2, lty=3, col="grey")
# Fc points to emphasize line
points(Fs_A1 ~ D_over_w
       , col= alpha(colorlist[2], point_transparency)
       , pch = 19
)
# N points to emphasize line
points(Ns_A1 ~ D_over_w
       , col= alpha(colorlist[1], point_transparency)
       , pch = 19
)
# Letter label for this panel
text(1-offsetx, Nlim*offsety, "a", cex=2, col = "black")
mtext(text = "Model A.1", 
      side = 3, line = 0, las=1, cex=1, xpd=TRUE)
# Text in local margin for this row's y-axis label
mtext(text = "Number of bees", 
      side = 2, line = 3, las=3, cex=0.8)
legend("topleft"
       , legend = c("N (total bees)", "N-F (bees on w)", "F (bees on c)", "Nc (bees on c at w=0")
       , col = c(colorlist[1], colorlist[5], colorlist[2], colorlist[3])
       , pch = 19
       , cex = 0.9
       , pt.cex = 2
       , bg="transparent"
       , bty = 'n'
       
)
# ROW 1: Number of bees: N, Fc, Fw=N-Fc, Nc
# column 2: Model A2
plot(NULL
     , ylim= c(0, Nlim)
     , xlim= c(0, 1)
     , xaxt= 'n'
     , yaxt= 'n'
)
# Fc filled area
polygon(c(D_over_w, rev(D_over_w)), c(c(1:resolution)*0, rev(Fs_A2))
        , col = alpha(colorlist[2], point_transparency)
        , lty = 0
)
# Fw = area between Fc and N
polygon(c(D_over_w, rev(D_over_w)), c(Fs_A2, rev(Ns_A2))
        , col = colorlist[5]
        , lty = 0
)
# Letter label for this panel
text(1-offsetx, Nlim*offsety, "b", cex=2, col = "black")
mtext(text = "Model A.2", 
      side = 3, line = 0, las=1, cex=1, xpd=TRUE)
# horizontal line at Nc
abline(h=N_A2(c_fixed, 0), lty=2, lwd=3, col=colorlist[3])
# Cross of lines at w*
abline(h=max(Fs_A2), lwd=2, lty=3, col="grey")
abline(v=D_over_w[which(Fs_A2==max(Fs_A2))], lwd=2, lty=3, col="grey")
# Fc points to emphasize line
points(Fs_A2 ~ D_over_w
       , col= alpha(colorlist[2], point_transparency)
       , pch = 19
)
# N points to emphasize line
points(Ns_A2 ~ D_over_w
       , col= alpha(colorlist[1], point_transparency)
       , pch = 19
)

# ROW 1: Number of bees: N, Fc, Fw=N-Fc, Nc
# column 3: Model A3
plot(NULL
     , ylim= c(0, Nlim)
     , xlim= c(0, 1)
     , xaxt= 'n'
     , yaxt= 'n'
)
# Fc filled area
polygon(c(D_over_w, rev(D_over_w)), c(c(1:resolution)*0, rev(Fs_A3))
        , col = alpha(colorlist[2], point_transparency)
        , lty = 0
)
# Fw = area between Fc and N
polygon(c(D_over_w, rev(D_over_w)), c(Fs_A3, rev(Ns_A3))
        , col = colorlist[5]
        , lty = 0
)
# Letter label for this panel
text(1-offsetx, Nlim*offsety, "c", cex=2, col = "black")
mtext(text = "Model A.3", 
      side = 3, line = 0, las=1, cex=1, xpd=TRUE)
# horizontal line at Nc
abline(h=N_A3(c_fixed, 0), lty=2, lwd=3, col=colorlist[3])
# Cross of lines at w*
abline(h=max(Fs_A3), lwd=2, lty=3, col="grey")
abline(v=D_over_w[which(Fs_A3==max(Fs_A3))], lwd=2, lty=3, col="grey")
# Fc points to emphasize line
points(Fs_A3 ~ D_over_w
       , col= alpha(colorlist[2], point_transparency)
       , pch = 19
)
# N points to emphasize line
points(Ns_A3 ~ D_over_w
       , col= alpha(colorlist[1], point_transparency)
       , pch = 19
)


# 2nd ROW - Differences caused by other species -----------------
# ROW 2: Fc-Nc (gain for c from w), Fw-N(0,w) (gain for w from c)
# column 1: Model A1
plot(NULL
     , ylim= c(-Nlim, Nlim)
     , xlim= c(0, 1)
     , xaxt= 'n'
)
# horizontal line at 0
abline(h=0, lwd=2, lty=2, col="grey")
# vertical lines at w*
abline(v=D_over_w[which(Fs_A1==max(Fs_A1))], lwd=2, lty=3, col="grey")
# Fc-Nc (gain for c from w)
points(Fs_A1-N_A1(c_fixed,0) ~ D_over_w
       , col= alpha(colorlist[2], point_transparency)
       , pch = 19
)
# Fw-N(0,w) (gain for w from c)
points((Ns_A1-Fs_A1)-N_A1(0,w_list) ~ D_over_w
       , col= alpha(colorlist[5], point_transparency)
       , pch = 19
)
# Letter label for this panel
text(1*offsetx, Nlim *(offsety*2-1), "c", cex=2, col = "black")
# Text in local margin for this row's y-axis label
mtext(text = "Bees gained or lost", 
      side = 2, line = 3, las=3, cex=0.8)

# ROW 2: Fc-Nc (gain for c from w), Fw-N(0,w) (gain for w from c)
# column 2: Model A2
plot(NULL
     , ylim= c(-Nlim, Nlim)
     , xlim= c(0, 1)
     , xaxt= 'n'
     , yaxt= 'n'
)
# horizontal line at 0
abline(h=0, lwd=2, lty=2, col="grey")
# vertical lines at w*
abline(v=D_over_w[which(Fs_A2==max(Fs_A2))], lwd=2, lty=3, col="grey")
# Fc-Nc (gain for c from w)
points(Fs_A2-N_A2(c_fixed,0) ~ D_over_w
       , col= alpha(colorlist[2], point_transparency)
       , pch = 19
)
# Fw-N(0,w) (gain for w from c)
points((Ns_A2-Fs_A2)-N_A2(0,w_list) ~ D_over_w
       , col= alpha(colorlist[5], point_transparency)
       , pch = 19
)
# Letter label for this panel
text(1*offsetx, Nlim *(offsety*2-1), "d", cex=2, col = "black")

# ROW 2: Fc-Nc (gain for c from w), Fw-N(0,w) (gain for w from c)
# column 3: Model A3
plot(NULL
     , ylim= c(-Nlim, Nlim)
     , xlim= c(0, 1)
     , xaxt= 'n'
     , yaxt= 'n'
)
# horizontal line at 0
abline(h=0, lwd=2, lty=2, col="grey")
# vertical lines at w*
abline(v=D_over_w[which(Fs_A3==max(Fs_A3))], lwd=2, lty=3, col="grey")
# Fc-Nc (gain for c from w)
points(Fs_A3-N_A3(c_fixed,0) ~ D_over_w
       , col= alpha(colorlist[2], point_transparency)
       , pch = 19
)
# Fw-N(0,w) (gain for w from c)
points((Ns_A3-Fs_A3)-N_A3(0,w_list) ~ D_over_w
       , col= alpha(colorlist[5], point_transparency)
       , pch = 19
)
# Letter label for this panel
text(1*offsetx, Nlim *(offsety*2-1), "e", cex=2, col = "black")
legend("topright"
       , legend = c("Fw-N(0,w) \n(bees gained by w because of c)", "F-Nc \n(bees gained by c because of w)")
       , col = c(colorlist[5], colorlist[2])
       , pch = 19
       , cex = 0.8
       , pt.cex = 2
       , bg="transparent"
       , bty = "n"
       , y.intersp = 1.5
)

# 3rd ROW - Bees per amount of c and w ---------------
# ROW 3: Fc/c and Fw/w
# column 1: Model A.1
plot(NULL
     , ylim= c(0, F_clim)
     , xlim= c(0, 1)
)
# vertical lines at w*
abline(v=D_over_w[which(Fs_A1==max(Fs_A1))], lwd=2, lty=3, col="grey")
# Fc/c 
points(Fs_A1/c_fixed ~ D_over_w
       , col= alpha(colorlist[2], point_transparency)
       , pch = 19
)
# Fw/w
points((Ns_A1-Fs_A1)/w_list ~ D_over_w
       , col= alpha(colorlist[5], point_transparency)
       , pch = 19
)
# Letter label for this panel
text(1*offsetx, F_clim*offsety, "g", cex=2, col = "black")
# Text in local margin for this row's y-axis label
mtext(text = "Bees per 'amount' of flowers", 
      side = 2, line = 3, las=3, cex=0.8)
#legend("topright"
#       , legend = c("Fw/w", "Fc/c")
#       , col = c(colorlist[5], colorlist[2])
#       , pch = 19
#       , cex = 0.9
#       , pt.cex = 2
#       #, bg="transparent"
#)

# ROW 3: Fc/c and Fw/w
# column 2: Model A2
plot(NULL
     , ylim= c(0, F_clim)
     , xlim= c(0, 1)
     , yaxt= 'n'
)
# vertical lines at w*
abline(v=D_over_w[which(Fs_A2==max(Fs_A2))], lwd=2, lty=3, col="grey")
# Fc/c 
points(Fs_A2/c_fixed ~ D_over_w
       , col= alpha(colorlist[2], point_transparency)
       , pch = 19
)
# Fw/w
points((Ns_A2-Fs_A2)/w_list ~ D_over_w
       , col= alpha(colorlist[5], point_transparency)
       , pch = 19
)
# Letter label for this panel
text(1*offsetx, F_clim*offsety, "h", cex=2, col = "black")

# ROW 3: Fc/c and Fw/w
# column 3: Model A3
plot(NULL
     , ylim= c(0, F_clim)
     , xlim= c(0, 1)
     , yaxt= 'n'
)
# vertical lines at w*
abline(v=D_over_w[which(Fs_A3==max(Fs_A3))], lwd=2, lty=3, col="grey")
# Fc/c 
points(Fs_A3/c_fixed ~ D_over_w
       , col= alpha(colorlist[2], point_transparency)
       , pch = 19
)
# Fw/w
points((Ns_A3-Fs_A3)/w_list ~ D_over_w
       , col= alpha(colorlist[5], point_transparency)
       , pch = 19
)
# Letter label for this panel
text(1*offsetx, F_clim*offsety, "i", cex=2, col = "black")

# X-Axis label ----------------------------
# Text in outer margin for overarching x-axis label
mtext(text = "Fraction of pollinators in patch on crop: D = c/(c+w)", 
      side = 1, line = 3, las=1, cex=0.8, outer=TRUE)

### End Fig S1 ----------------------

### FIG S2 PARAMETERS ON D, and w* 3x3 ------------------------
# x-axis D, y-axis N, 3 columns = 3 models
# 4 rows: c, Aw, aw, ac in colors
# Figure format ---------------------
pointsize <- 0.6
Nlim <- 20 # Max for N graph
F_clim <- Nlim / 10 # Max for F/c graph
wlim <- 20 # for plotting
clim <- 50 # Max c for graph
offsetx <- 0.08
offsety <- 0.9
point_transparency <- 0.5
opt_w_lines <- 3
F_lines <- 1
N_lines <- 2
leg_x <- -0.7
leg_y <- 0.9

# Make 3x4 panels
par(mfcol=c(4,3))
# Various other format adjustments (margins in order bottom, left, top, right)
par(oma = c(4,12,1,0), mar=c(0.5, 0.5, 0.5, 0.5), mgp=c(3, 1, 0), las=1, xpd=FALSE)
# Parameters --------------------
default_c <- c_fixed
default_AAw <- AAw
default_aw <- aw
default_ac <- ac
linewidth <- 4
linetransp <- 0.8

# 1st COLUMN Model A1 --------------------
# PANEL A: A.1, varying over c -------------------
# Initialize 2-d arrays
# Each is a list over w, but with rows being parameter state
N_figS2 <- array(NaN, c(no_states, resolution))
F_figS2 <- array(NaN, c(no_states, resolution))
Nc_figS2 <- array(NaN, c(no_states, resolution))
w_opt_figS2 <- array(NaN, c(no_states))

for (i in 1:no_states) { # over parameter values
  c_fixed <- c_states[i]
  D_over_w <- D(c_fixed, w_list)
  N_figS2[i,] <- N_A1(c_fixed,w_list)
  F_figS2[i,] <- F(N_figS2[i,], D_over_w)
  Nc_figS2[i,] <- N_A1(c_fixed,rep(0, resolution))
  w_opt_figS2[i] <- w_list[which(F_figS2[i,]==max(F_figS2[i,]))]
}
c_fixed <- default_c
D_over_w <- D(c_fixed, w_list)

plot(NULL
     , ylim= c(0, Nlim)
     , xlim= c(0, 1)
     , xaxt= 'n'
)
for (i in 1:no_states) { # over parameter values
  lines(N_figS2[i,]~D(c_states[i], w_list)
        , col = N_colors[i]
        , lwd = linewidth
        , lty = N_lines
  )
  lines(F_figS2[i,]~D(c_states[i], w_list)
        , col = F_colors[i]
        , lwd = linewidth
        , lty = F_lines
  )
  abline(v=D(c_states[i], w_opt_figS2[i]), lty = opt_w_lines, lwd=2, col=F_colors[i])
}

text(D(c_states[1], w_opt_figS2[1])-0.015, Nlim*0.7
       , paste("w*=", round(w_opt_figS2[1], 1))
       , cex=1, col="black", pos=2, xpd=NA
     )
text(D(c_states[no_states], w_opt_figS2[no_states])+0.015, Nlim*0.7
     , paste("w*=\n", round(w_opt_figS2[no_states], 1))
     , cex=1, col="black", pos=4, xpd=NA
)
# Letter label for this panel
text(offsetx, Nlim*offsety, "a", cex=2, col = "black")
mtext(text = "Model A.1", 
      side = 3, line = 0, las=1, cex=1, xpd=TRUE)
# Legends & Row titles in left margin -----------------------------
text(leg_x, Nlim, "Varying c \n(amount of crop)"
     , cex=1, col="black", pos=4, xpd=NA)
legend(x=leg_x
       , y=Nlim*leg_y
       , xpd = NA
       , legend = paste("c=", c_states)
       , text.col = "transparent"
       , col = N_colors
       , pch = 19
       , cex = 0.9
       , pt.cex = 2
       , bty = 'n'
       , bg= "transparent"
)
legend(x=leg_x+0.1
       , y=Nlim*leg_y
       , xpd = NA
       , legend = paste("c=", c_states)
       , col = F_colors
       , pch = 19
       , cex = 0.9
       , pt.cex = 2
       , bty = 'n'
       , bg= "transparent"
)
# PANEL D: A.1, varying over AAw -------------------

# Initialize 2-d arrays
# Each is a list over w, but with rows being parameter state
N_figS2 <- array(NaN, c(no_states, resolution))
F_figS2 <- array(NaN, c(no_states, resolution))
Nc_figS2 <- array(NaN, c(no_states, resolution))
w_opt_figS2 <- array(NaN, c(no_states))

for (i in 1:no_states) { # over parameter values
  AAw <- AAw_list[i]
  N_figS2[i,] <- N_A1(c_fixed,w_list)
  F_figS2[i,] <- F(N_figS2[i,], D_over_w)
  Nc_figS2[i,] <- N_A1(c_fixed,rep(0, resolution))
  w_opt_figS2[i] <- w_list[which(F_figS2[i,]==max(F_figS2[i,]))]
}
AAw <- default_AAw

plot(NULL
     , ylim= c(0, Nlim)
     , xlim= c(0, 1)
     , xaxt= 'n'
)
for (i in 1:no_states) { # over parameter values
  lines(N_figS2[i,]~D_over_w
        , col = N_colors[i]
        , lwd = linewidth
        , lty = N_lines
  )
  lines(F_figS2[i,]~D_over_w
        , col = F_colors[i]
        , lwd = linewidth
        , lty = F_lines
  )
  abline(v=D(c_fixed, w_opt_figS2[i]), lty = opt_w_lines, lwd=2, col=F_colors[i])
}

# Letter label for this panel
text(offsetx, Nlim*offsety, "d", cex=2, col = "black")
text(D(c_fixed, w_opt_figS2[no_states])-0.015, Nlim*0.8
     , paste("w*=", round(w_opt_figS2[no_states], 1))
     , cex=1, col="black", pos=2, xpd=NA
)
# Legends & Row titles in left margin -----------------------------
text(leg_x, Nlim, "Varying Aw \n(attractivenes of w)"
     , cex=1, col="black", pos=4, xpd=NA)
legend(x=leg_x
       , y=Nlim*leg_y
       , xpd = NA
       , legend = paste("c=", c_states)
       , text.col = "transparent"
       , col = N_colors
       , pch = 19
       , cex = 0.9
       , pt.cex = 2
       , bty = 'n'
       , bg= "transparent"
)
legend(x=leg_x+0.1
       , y=Nlim*leg_y
       , xpd = NA
       , legend = paste("Aw=", AAw_list)
       , col = F_colors
       , pch = 19
       , cex = 0.9
       , pt.cex = 2
       , bty = 'n'
       , bg= "transparent"
)
# PANEL G: A.1, varying over aw -------------------

# Initialize 2-d arrays
# Each is a list over w, but with rows being parameter state
N_figS2 <- array(NaN, c(no_states, resolution))
F_figS2 <- array(NaN, c(no_states, resolution))
Nc_figS2 <- array(NaN, c(no_states, resolution))
w_opt_figS2 <- array(NaN, c(no_states))

for (i in 1:no_states) { # over parameter values
  aw <- aw_list[i]
  N_figS2[i,] <- N_A1(c_fixed,w_list)
  F_figS2[i,] <- F(N_figS2[i,], D_over_w)
  Nc_figS2[i,] <- N_A1(c_fixed,rep(0, resolution))
  w_opt_figS2[i] <- w_list[which(F_figS2[i,]==max(F_figS2[i,]))]
}
aw <- default_aw

plot(NULL
     , ylim= c(0, Nlim)
     , xlim= c(0, 1)
     , xaxt= 'n'
)
for (i in 1:no_states) { # over parameter values
  lines(N_figS2[i,]~D_over_w
        , col = N_colors[i]
        , lwd = linewidth
        , lty = N_lines
  )
  lines(F_figS2[i,]~D_over_w
        , col = F_colors[i]
        , lwd = linewidth
        , lty = F_lines
  )
  abline(v=D(c_fixed, w_opt_figS2[i]), lty = opt_w_lines, lwd=2, col=F_colors[i])
}
# Letter label for this panel
text(offsetx, Nlim*offsety, "g", cex=2, col = "black")
text(D(c_fixed, w_opt_figS2[no_states])-0.015, Nlim*0.7
     , paste("w*=", round(w_opt_figS2[no_states], 1))
     , cex=1, col="black", pos=2, xpd=NA
)
# Legends & Row titles in left margin -----------------------------
text(leg_x, Nlim, "Varying aw \n(shape of N(w))"
     , cex=1, col="black", pos=4, xpd=NA)
legend(x=leg_x
       , y=Nlim*leg_y
       , xpd = NA
       , legend = paste("c=", c_states)
       , text.col = "transparent"
       , col = N_colors
       , pch = 19
       , cex = 0.9
       , pt.cex = 2
       , bty = 'n'
       , bg= "transparent"
)
legend(x=leg_x+0.1
       , y=Nlim*leg_y
       , xpd = NA
       , legend = paste("aw=", aw_list)
       , col = F_colors
       , pch = 19
       , cex = 0.9
       , pt.cex = 2
       , bty = 'n'
       , bg= "transparent"
)
# PANEL J: A.1, varying over ac -------------------

# Initialize 2-d arrays
# Each is a list over w, but with rows being parameter state
N_figS2 <- array(NaN, c(no_states, resolution))
F_figS2 <- array(NaN, c(no_states, resolution))
Nc_figS2 <- array(NaN, c(no_states, resolution))
w_opt_figS2 <- array(NaN, c(no_states))

for (i in 1:no_states) { # over parameter values
  ac <- ac_list[i]
  N_figS2[i,] <- N_A1(c_fixed,w_list)
  F_figS2[i,] <- F(N_figS2[i,], D_over_w)
  Nc_figS2[i,] <- N_A1(c_fixed,rep(0, resolution))
  w_opt_figS2[i] <- w_list[which(F_figS2[i,]==max(F_figS2[i,]))]
}
ac <- default_ac

plot(NULL
     , ylim= c(0, Nlim)
     , xlim= c(0, 1)
)
for (i in 1:no_states) { # over parameter values
  lines(N_figS2[i,]~D_over_w
        , col = N_colors[i]
        , lwd = linewidth
        , lty = N_lines
  )
  lines(F_figS2[i,]~D_over_w
        , col = F_colors[i]
        , lwd = linewidth
        , lty = F_lines
  )
  abline(v=D(c_fixed, w_opt_figS2[i]), lty = opt_w_lines, lwd=2, col=F_colors[i])
}
text(D(c_fixed, w_opt_figS2[no_states])-0.015, Nlim*0.7
     , paste("w*=", round(w_opt_figS2[no_states], 1))
     , cex=1, col="black", pos=2, xpd=NA
)
# Letter label for this panel
text(offsetx, Nlim*offsety, "j", cex=2, col = "black")
# Legends & Row titles in left margin -----------------------------
text(leg_x, Nlim, "Varying ac \n(shape of N(c))"
     , cex=1, col="black", pos=4, xpd=NA)
legend(x=leg_x
       , y=Nlim*leg_y
       , xpd = NA
       , legend = paste("c=", c_states)
       , text.col = "transparent"
       , col = N_colors
       , pch = 19
       , cex = 0.9
       , pt.cex = 2
       , bty = 'n'
       , bg= "transparent"
)
legend(x=leg_x+0.1
       , y=Nlim*leg_y
       , xpd = NA
       , legend = paste("ac=", ac_list)
       , col = F_colors
       , pch = 19
       , cex = 0.9
       , pt.cex = 2
       , bty = 'n'
       , bg= "transparent"
)
legend(x=-0.85
       , y=Nlim*0.16
       , xpd = NA
       , legend = c("N: Total bees \nin patch", "F: bees on c")
       , col = c(N_colors[3],F_colors[3])
       , lwd = 2
       , lty = c(N_lines, F_lines)
       , cex = 1.2
       , bty = 'n'
       , bg= "transparent"
)
# 2nd COLUMN Model A2 --------------------
# PANEL B: A.2, varying over c -------------------
# Initialize 2-d arrays
# Each is a list over w, but with rows being parameter state
N_figS2 <- array(NaN, c(no_states, resolution))
F_figS2 <- array(NaN, c(no_states, resolution))
Nc_figS2 <- array(NaN, c(no_states, resolution))
w_opt_figS2 <- array(NaN, c(no_states))

for (i in 1:no_states) { # over parameter values
  c_fixed <- c_states[i]
  D_over_w <- D(c_fixed, w_list)
  N_figS2[i,] <- N_A2(c_fixed,w_list)
  F_figS2[i,] <- F(N_figS2[i,], D_over_w)
  Nc_figS2[i,] <- N_A2(c_fixed,rep(0, resolution))
  w_opt_figS2[i] <- w_list[which(F_figS2[i,]==max(F_figS2[i,]))]
}
c_fixed <- default_c
D_over_w <- D(c_fixed, w_list)

plot(NULL
     , ylim= c(0, Nlim)
     , xlim= c(0, 1)
     , xaxt= 'n'
     , yaxt= 'n'
)
for (i in 1:no_states) { # over parameter values
  lines(N_figS2[i,]~D(c_states[i], w_list)
        , col = N_colors[i]
        , lwd = linewidth
        , lty = N_lines
  )
  lines(F_figS2[i,]~D(c_states[i], w_list)
        , col = F_colors[i]
        , lwd = linewidth
        , lty = F_lines
  )
  abline(v=D(c_states[i], w_opt_figS2[i]), lty = opt_w_lines, lwd=2, col=F_colors[i])
}
text(D(c_states[1], w_opt_figS2[1])-0.3, Nlim*0.1
     , paste("w*=", round(w_opt_figS2[1], 1))
     , cex=1, col=F_colors[1], pos=4, xpd=NA
)
text(D(c_states[no_states], w_opt_figS2[no_states])+0.06, Nlim*0.1
     , paste("w*=", round(w_opt_figS2[no_states], 1))
     , cex=1, col=F_colors[no_states], pos=4, xpd=NA
)
# Letter label for this panel
text(offsetx, Nlim*offsety, "b", cex=2, col = "black")
mtext(text = "Model A.2", 
      side = 3, line = 0, las=1, cex=1, xpd=TRUE)

# PANEL E: A.2, varying over AAw -------------------

# Initialize 2-d arrays
# Each is a list over w, but with rows being parameter state
N_figS2 <- array(NaN, c(no_states, resolution))
F_figS2 <- array(NaN, c(no_states, resolution))
Nc_figS2 <- array(NaN, c(no_states, resolution))
w_opt_figS2 <- array(NaN, c(no_states))

for (i in 1:no_states) { # over parameter values
  AAw <- AAw_list[i]
  N_figS2[i,] <- N_A2(c_fixed,w_list)
  F_figS2[i,] <- F(N_figS2[i,], D_over_w)
  Nc_figS2[i,] <- N_A2(c_fixed,rep(0, resolution))
  w_opt_figS2[i] <- w_list[which(F_figS2[i,]==max(F_figS2[i,]))]
}
AAw <- default_AAw

plot(NULL
     , ylim= c(0, Nlim)
     , xlim= c(0, 1)
     , xaxt= 'n'
     , yaxt= 'n'
)
for (i in 1:no_states) { # over parameter values
  lines(N_figS2[i,]~D_over_w
        , col = N_colors[i]
        , lwd = linewidth
        , lty = N_lines
  )
  lines(F_figS2[i,]~D_over_w
        , col = F_colors[i]
        , lwd = linewidth
        , lty = F_lines
  )
  abline(v=D(c_fixed, w_opt_figS2[i]), lty = opt_w_lines, lwd=2, col=F_colors[i])
}
text(D(c_fixed, w_opt_figS2[no_states])-0.015, Nlim*0.8
     , paste("w*=", round(w_opt_figS2[no_states], 1))
     , cex=1, col="black", pos=2, xpd=NA
)
# Letter label for this panel
text(offsetx, Nlim*offsety, "e", cex=2, col = "black")

# PANEL H: A.2, varying over aw -------------------

# Initialize 2-d arrays
# Each is a list over w, but with rows being parameter state
N_figS2 <- array(NaN, c(no_states, resolution))
F_figS2 <- array(NaN, c(no_states, resolution))
Nc_figS2 <- array(NaN, c(no_states, resolution))
w_opt_figS2 <- array(NaN, c(no_states))

for (i in 1:no_states) { # over parameter values
  aw <- aw_list[i]
  N_figS2[i,] <- N_A2(c_fixed,w_list)
  F_figS2[i,] <- F(N_figS2[i,], D_over_w)
  Nc_figS2[i,] <- N_A2(c_fixed,rep(0, resolution))
  w_opt_figS2[i] <- w_list[which(F_figS2[i,]==max(F_figS2[i,]))]
}
aw <- default_aw

plot(NULL
     , ylim= c(0, Nlim)
     , xlim= c(0, 1)
     , yaxt= 'n'
     , xaxt= 'n'
)
for (i in 1:no_states) { # over parameter values
  lines(N_figS2[i,]~D_over_w
        , col = N_colors[i]
        , lwd = linewidth
        , lty = N_lines
  )
  lines(F_figS2[i,]~D_over_w
        , col = F_colors[i]
        , lwd = linewidth
        , lty = F_lines
  )
  abline(v=D(c_fixed, w_opt_figS2[i]), lty = opt_w_lines, lwd=2, col=F_colors[i])
}
# Letter label for this panel
text(offsetx, Nlim*offsety, "h", cex=2, col = "black")

# PANEL K: A.2, varying over ac -------------------

# Initialize 2-d arrays
# Each is a list over w, but with rows being parameter state
N_figS2 <- array(NaN, c(no_states, resolution))
F_figS2 <- array(NaN, c(no_states, resolution))
Nc_figS2 <- array(NaN, c(no_states, resolution))
w_opt_figS2 <- array(NaN, c(no_states))

for (i in 1:no_states) { # over parameter values
  ac <- ac_list[i]
  N_figS2[i,] <- N_A2(c_fixed,w_list)
  F_figS2[i,] <- F(N_figS2[i,], D_over_w)
  Nc_figS2[i,] <- N_A2(c_fixed,rep(0, resolution))
  w_opt_figS2[i] <- w_list[which(F_figS2[i,]==max(F_figS2[i,]))]
}
ac <- default_ac

plot(NULL
     , ylim= c(0, Nlim)
     , xlim= c(0, 1)
     , yaxt= 'n'
)
for (i in 1:no_states) { # over parameter values
  lines(N_figS2[i,]~D_over_w
        , col = N_colors[i]
        , lwd = linewidth
        , lty = N_lines
  )
  lines(F_figS2[i,]~D_over_w
        , col = F_colors[i]
        , lwd = linewidth
        , lty = F_lines
  )
  abline(v=D(c_fixed, w_opt_figS2[i]), lty = opt_w_lines, lwd=2, col=F_colors[i])
}
text(D(c_fixed, w_opt_figS2[1])-0.015, Nlim*0.2
     , paste("w*=", round(w_opt_figS2[1], 1))
     , cex=1, col="black", pos=2, xpd=NA
)
text(D(c_fixed, w_opt_figS2[no_states])+0.015, Nlim*0.2
     , paste("w*=", round(w_opt_figS2[no_states], 1))
     , cex=1, col="black", pos=4, xpd=NA
)
# Letter label for this panel
text(offsetx, Nlim*offsety, "k", cex=2, col = "black")

# 3rd COLUMN Model A3 --------------------
# PANEL C: A.3, varying over c -------------------
# Initialize 2-d arrays
# Each is a list over w, but with rows being parameter state
N_figS2 <- array(NaN, c(no_states, resolution))
F_figS2 <- array(NaN, c(no_states, resolution))
Nc_figS2 <- array(NaN, c(no_states, resolution))
w_opt_figS2 <- array(NaN, c(no_states))

for (i in 1:no_states) { # over parameter values
  c_fixed <- c_states[i]
  D_over_w <- D(c_fixed, w_list)
  N_figS2[i,] <- N_A3(c_fixed,w_list)
  F_figS2[i,] <- F(N_figS2[i,], D_over_w)
  Nc_figS2[i,] <- N_A3(c_fixed,rep(0, resolution))
  w_opt_figS2[i] <- w_list[which(F_figS2[i,]==max(F_figS2[i,]))]
}
c_fixed <- default_c
D_over_w <- D(c_fixed, w_list)

plot(NULL
     , ylim= c(0, Nlim)
     , xlim= c(0, 1)
     , xaxt= 'n'
     , yaxt= 'n'
)
for (i in 1:no_states) { # over parameter values
  lines(N_figS2[i,]~D(c_states[i], w_list)
        , col = N_colors[i]
        , lwd = linewidth
        , lty = N_lines
  )
  lines(F_figS2[i,]~D(c_states[i], w_list)
        , col = F_colors[i]
        , lwd = linewidth
        , lty = F_lines
  )
  abline(v=D(c_states[i], w_opt_figS2[i]), lty = opt_w_lines, lwd=2, col=F_colors[i])
}
text(D(c_states[1], w_opt_figS2[1])-0.02, Nlim*0.7
     , paste("w*=", round(w_opt_figS2[1], 1))
     , cex=1, col="black", pos=2, xpd=NA
)
text(D(c_states[no_states], w_opt_figS2[no_states])+0.02, Nlim*0.7
     , paste("w*=", round(w_opt_figS2[no_states], 1))
     , cex=1, col="black", pos=4, xpd=NA
)
# Letter label for this panel
text(offsetx, Nlim*offsety, "c", cex=2, col = "black")
mtext(text = "Model A.3", 
      side = 3, line = 0, las=1, cex=1, xpd=TRUE)

# PANEL F: A.3, varying over AAw -------------------

# Initialize 2-d arrays
# Each is a list over w, but with rows being parameter state
N_figS2 <- array(NaN, c(no_states, resolution))
F_figS2 <- array(NaN, c(no_states, resolution))
Nc_figS2 <- array(NaN, c(no_states, resolution))
w_opt_figS2 <- array(NaN, c(no_states))

for (i in 1:no_states) { # over parameter values
  AAw <- AAw_list[i]
  N_figS2[i,] <- N_A3(c_fixed,w_list)
  F_figS2[i,] <- F(N_figS2[i,], D_over_w)
  Nc_figS2[i,] <- N_A3(c_fixed,rep(0, resolution))
  w_opt_figS2[i] <- w_list[which(F_figS2[i,]==max(F_figS2[i,]))]
}
AAw <- default_AAw

plot(NULL
     , ylim= c(0, Nlim)
     , xlim= c(0, 1)
     , xaxt= 'n'
     , yaxt= 'n'
)
for (i in 1:no_states) { # over parameter values
  lines(N_figS2[i,]~D_over_w
        , col = N_colors[i]
        , lwd = linewidth
        , lty = N_lines
  )
  lines(F_figS2[i,]~D_over_w
        , col = F_colors[i]
        , lwd = linewidth
        , lty = F_lines
  )
  abline(v=D(c_fixed, w_opt_figS2[i]), lty = opt_w_lines, lwd=2, col=F_colors[i])
}
text(D(c_fixed, w_opt_figS2[no_states])-0.25, Nlim*0.8
     , paste("w*=", round(w_opt_figS2[no_states], 1))
     , cex=1, col="black", pos=4, xpd=NA
)
# Letter label for this panel
text(offsetx, Nlim*offsety, "f", cex=2, col = "black")

# PANEL I: A.3, varying over aw -------------------

# Initialize 2-d arrays
# Each is a list over w, but with rows being parameter state
N_figS2 <- array(NaN, c(no_states, resolution))
F_figS2 <- array(NaN, c(no_states, resolution))
Nc_figS2 <- array(NaN, c(no_states, resolution))
w_opt_figS2 <- array(NaN, c(no_states))

for (i in 1:no_states) { # over parameter values
  aw <- aw_list[i]
  N_figS2[i,] <- N_A3(c_fixed,w_list)
  F_figS2[i,] <- F(N_figS2[i,], D_over_w)
  Nc_figS2[i,] <- N_A3(c_fixed,rep(0, resolution))
  w_opt_figS2[i] <- w_list[which(F_figS2[i,]==max(F_figS2[i,]))]
}
aw <- default_aw

plot(NULL
     , ylim= c(0, Nlim)
     , xlim= c(0, 1)
     , yaxt= 'n'
     , xaxt= 'n'
)
for (i in 1:no_states) { # over parameter values
  lines(N_figS2[i,]~D_over_w
        , col = N_colors[i]
        , lwd = linewidth
        , lty = N_lines
  )
  lines(F_figS2[i,]~D_over_w
        , col = F_colors[i]
        , lwd = linewidth
        , lty = F_lines
  )
  abline(v=D(c_fixed, w_opt_figS2[i]), lty = opt_w_lines, lwd=2, col=F_colors[i])
}
text(D(c_fixed, w_opt_figS2[no_states])+0.015, Nlim*0.7
     , paste("w*=\n", round(w_opt_figS2[no_states], 1))
     , cex=1, col="black", pos=4, xpd=NA
)
text(D(c_fixed, w_opt_figS2[2])-0.015, Nlim*0.7
     , paste("w*=", round(w_opt_figS2[2], 1))
     , cex=1, col="black", pos=2, xpd=NA
)
# Letter label for this panel
text(offsetx, Nlim*offsety, "i", cex=2, col = "black")
# PANEL L: A.3, varying over ac -------------------

# Initialize 2-d arrays
# Each is a list over w, but with rows being parameter state
N_figS2 <- array(NaN, c(no_states, resolution))
F_figS2 <- array(NaN, c(no_states, resolution))
Nc_figS2 <- array(NaN, c(no_states, resolution))
w_opt_figS2 <- array(NaN, c(no_states))

for (i in 1:no_states) { # over parameter values
  ac <- ac_list[i]
  N_figS2[i,] <- N_A3(c_fixed,w_list)
  F_figS2[i,] <- F(N_figS2[i,], D_over_w)
  Nc_figS2[i,] <- N_A3(c_fixed,rep(0, resolution))
  w_opt_figS2[i] <- w_list[which(F_figS2[i,]==max(F_figS2[i,]))]
}
ac <- default_ac

plot(NULL
     , ylim= c(0, Nlim)
     , xlim= c(0, 1)
     , yaxt= 'n'
)
for (i in 1:no_states) { # over parameter values
  lines(N_figS2[i,]~D_over_w
        , col = N_colors[i]
        , lwd = linewidth
        , lty = N_lines
  )
  lines(F_figS2[i,]~D_over_w
        , col = F_colors[i]
        , lwd = linewidth
        , lty = F_lines
  )
  abline(v=D(c_fixed, w_opt_figS2[i]), lty = opt_w_lines, lwd=2, col=F_colors[i])
}
text(D(c_fixed, w_opt_figS2[no_states])+0.02, Nlim*0.7
     , paste("w*=", round(w_opt_figS2[no_states], 1))
     , cex=1, col="black", pos=4, xpd=NA
)
# Letter label for this panel
text(offsetx, Nlim*offsety, "l", cex=2, col = "black")

# X-Axis label ----------------------------
# Text in outer margin for overarching x-axis label
mtext(text = "Fraction of pollinators in patch on crop: D = c/(c+w)", 
      side = 1, line = 2, las=1, cex=0.8, outer=TRUE)
# Y-Axis label ----------------------------
# Text in outer margin for overarching y-axis label
#mtext(text = "Total bees in patch N", 
#      side = 2, line = 2, las=3, cex=0.8, outer=TRUE)
### END FIG S2 --------------------------




### MODEL B ------------------
### GENERAL FUNCTIONS -----------------
c_B <- function(p) {
  return(1-p)
}
w_B <- function(p, b) {
  return(b*p)
}
crop_area <- function(p, total_area) {
  return(total_area*(1-p))
}
N_eq14 <- function(p, b, total_area, AAw, AAc, aw, ac) {
  return(N_eq14_c(p, b, total_area, AAw, AAc, aw, ac) + N_eq14_w(p, b, total_area, AAw, AAc, aw, ac))
}
N_eq14_c <- function(p, b, total_area, AAw, AAc, aw, ac) {
  return(AAc*(1-p) / (ac+1-p))
}
N_eq14_w <- function(p, b, total_area, AAw, AAc, aw, ac) {
  return(AAw*b*p / (aw+b*p))
}

D_B <- function(p, b) {
  return(c_B(p)/(c_B(p) + w_B(p, b)))
}
yield_B <- function(p, b, total_area, AAw, AAc, aw, ac) {
  return(N_eq14(p, b, total_area, AAw, AAc, aw, ac)*D_B(p, b)*crop_area(p, total_area))
}
# Function for graph later:

hatched_barplot_horizontal <- function(heights, names = NULL, hatch_spacing = 0.2, bar_height = 0.8, xlim = NULL, ...) {
  if (is.null(xlim)) xlim <- c(0, max(heights) * 1.1)
  bar_locs <- barplot(heights
                      , names.arg = names
                      , horiz = TRUE
                      , col = N_colors
                      , yaxt = 'n'
                      , xaxt = 'n'
                      , xpd = FALSE
                      , border = "black", xlim = c(0, xbars), ...)
  
  for (i in seq_along(bar_locs)) {
    y_center <- bar_locs[i]
    y_bottom <- y_center - bar_height / 2
    y_top <- y_center + bar_height / 2
    x_right <- heights[i]
    
    for (y in seq(y_bottom, y_top, by = hatch_spacing)) {
      lines(x = c(0, x_right), y = c(y, y), col = "white", lwd = 2)
    }
  }
  
  invisible(bar_locs)
}

# Parameters overall ----------------------
b <- 1
total_area <- 1
aw <- 0.4 
AAw <- 1
AAc <- 1
ac <- 0.4

# Over the range of p now repeating everything 
# with parameter variations; for each parameter, a new matrix. 
resolution_B <- 1000
p_vector <- seq(0, 1, length.out=resolution_B)
b_vector <- seq(1,100, length.out=no_states)
expon <- seq(-1, 1, length.out=no_states)
b_vector <- round(10^expon, 2)

### FIG S3 PARAMETERS ON p* 2x2 ------------------------
# x-axis D, y-axis y, 3 columns = 3 models
# 3 rows: Aw, aw, ac in colors, M and qa, qc for model A.3
# Figure format ---------------------
pointsize <- 0.6
offsetx <- 0.12
offsety <- 0.9
point_transparency <- 0.5
linewidth <- 4
linetransp <- 0.8
p_max <- 0.5
Ylim_B <- 2 # Max yield for graph
Y_clim_B <- 2 # Max for yield/c graph
leg_x <- -0.7
xbars <- 1

# Make 4x3 panels
par(mfrow=c(4,2))
# Various other format adjustments (margins in order bottom, left, top, right)
par(oma = c(5,12,0.5,5), mar=c(0.5, 0.5, 0.5, 0.5), mgp=c(3, 1, 0), las=1)
# Parameters --------------------
default_b <- b
default_AAw <- AAw
default_aw <- aw
default_ac <- ac

# ROW 1: B, varying over b -------------------
N_figS3 <- array(NaN, c(no_states, resolution_B))
Y_figS3 <- array(NaN, c(no_states, resolution_B))
p_opt_figS3 <- array(NaN, c(no_states))

for (i in 1:no_states) { # over parameter values
  b <- b_vector[i]
  D_over_p <- D_B(p_vector, b)
  N_figS3[i,] <- N_eq14(p_vector, b, total_area, AAw, AAc, aw, ac)
  Y_figS3[i,] <- yield_B(p_vector, b, total_area, AAw, AAc, aw, ac)
  p_opt_figS3[i] <- p_vector[which( Y_figS3[i,]==max(Y_figS3[i,]))]
}
b <- default_b
D_over_p <- D_B(p_vector, b)

plot(NULL
     , ylim= c(0, Ylim_B)
     , xlim= c(0, 1)
     , xaxt= 'n'
)
for (i in 1:no_states) { # over parameter values
  lines(Y_figS3[i,]~D_B(p_vector, b_vector[i])
        , col = F_colors[i]
        , lwd = linewidth
        , lty = F_lines
  )
  lines(N_figS3[i,]~D_B(p_vector, b_vector[i])
        , col = N_colors[i]
        , lwd = linewidth
        , lty = N_lines
  )
  abline(v=D_B(p_opt_figS3[i], b_vector[i]), lty = opt_w_lines, lwd=2, col=F_colors[i])
}
# Letter label for this panel
text(offsetx, Ylim_B*offsety, "a", cex=2, col = "black")
# Legends & Row titles in left margin -----------------------------
text(leg_x, Ylim_B, "Varying b \n(quality of w)"
     , cex=1, col="black", pos=4, xpd=NA)
legend(x=leg_x
       , y=Ylim_B*leg_y
       , xpd = NA
       , legend = paste("c=", c_states)
       , text.col = "transparent"
       , col = N_colors
       , pch = 19
       , cex = 0.9
       , pt.cex = 2
       , bty = 'n'
       , bg= "transparent"
)
legend(x=leg_x+0.1
       , y=Ylim_B*leg_y
       , xpd = NA
       , legend = paste("b=", b_vector)
       , col = F_colors
       , pch = 19
       , cex = 0.9
       , pt.cex = 2
       , bty = 'n'
       , bg= "transparent"
)
# 2nd COLUMN ROW 1 -------------------------------
hatched_barplot_horizontal(rep(xbars, no_states))
#barplot(rep(xbars, no_states)~b_vector
#        , col = N_colors
#        , xlim = c(0, xbars)
#        , horiz = T
#        , yaxt = 'n'
#        , xaxt = 'n'
#        , xpd = FALSE
#)
barplot((1-p_opt_figS3)~b_vector
        , col = F_colors
        , xlim = c(0, xbars)
        , horiz = T
        , yaxt = 'n'
        , xaxt = 'n'
        , xpd = FALSE
        , add = TRUE
        , lwd = 2
        )
box(lwd = 1)
# Letter label for this panel
text(xbars*(1-offsetx), 5.5, "b", cex=2, col = "black")
p_opt_figS3 <- round(p_opt_figS3, 2)
p_star_vals <- paste("p*=", p_opt_figS3[5], "\n\n"
                     , "p*=", p_opt_figS3[4], "\n\n"
                     , "p*=", p_opt_figS3[3], "\n\n"
                     , "p*=", p_opt_figS3[2], "\n\n"
                     , "p*=", p_opt_figS3[1], sep = ""
                     )
mtext(text = p_star_vals, 
      side = 4, line = 1, las=1, cex=0.8, outer=FALSE)

# ROW 2: B, varying over AAw -------------------
N_figS3 <- array(NaN, c(no_states, resolution_B))
Y_figS3 <- array(NaN, c(no_states, resolution_B))
p_opt_figS3 <- array(NaN, c(no_states))
D_over_p <- D_B(p_vector, b)

for (i in 1:no_states) { # over parameter values
  AAw <- AAw_list[i]
  N_figS3[i,] <- N_eq14(p_vector, b, total_area, AAw, AAc, aw, ac)
  Y_figS3[i,] <- yield_B(p_vector, b, total_area, AAw, AAc, aw, ac)
  p_opt_figS3[i] <- p_vector[which(Y_figS3[i,]==max(Y_figS3[i,]))]
}
AAw <- default_AAw

plot(NULL
     , ylim= c(0, Ylim_B)
     , xlim= c(0, 1)
     , xaxt='n'
)
for (i in 1:no_states) { # over parameter values
  lines(Y_figS3[i,]~D_B(p_vector, b)
        , col = F_colors[i]
        , lwd = linewidth
        , lty = F_lines
  )
  lines(N_figS3[i,]~D_B(p_vector, b)
        , col = N_colors[i]
        , lwd = linewidth
        , lty = N_lines
  )
  abline(v=D_B(p_opt_figS3[i], b), lty = opt_w_lines, lwd=2, col=F_colors[i])
}
# Letter label for this panel
text(offsetx, Ylim_B*offsety, "c", cex=2, col = "black")
# Legends & Row titles in left margin -----------------------------
text(leg_x, Ylim_B, "Varying Aw \n(attractivenes of w)"
     , cex=1, col="black", pos=4, xpd=NA)
legend(x=leg_x
       , y=Ylim_B*leg_y
       , xpd = NA
       , legend = paste("c=", c_states)
       , text.col = "transparent"
       , col = N_colors
       , pch = 19
       , cex = 0.9
       , pt.cex = 2
       , bty = 'n'
       , bg= "transparent"
)
legend(x=leg_x+0.1
       , y=Ylim_B*leg_y
       , xpd = NA
       , legend = paste("Aw=", AAw_list)
       , col = F_colors
       , pch = 19
       , cex = 0.9
       , pt.cex = 2
       , bty = 'n'
       , bg= "transparent"
)
# 2nd COLUMN ROW 2 -------------------------------
hatched_barplot_horizontal(rep(xbars, no_states))
barplot((1-p_opt_figS3)~b_vector
        , col = F_colors
        , xlim = c(0, xbars)
        , horiz = T
        , yaxt = 'n'
        , xaxt = 'n'
        , xpd = FALSE
        , add = TRUE
        , lwd = 2
)
box(lwd = 1)
# Letter label for this panel
text(xbars*(1-offsetx), 5.5, "d", cex=2, col = "black")
p_opt_figS3 <- round(p_opt_figS3, 2)
p_star_vals <- paste("p*=", p_opt_figS3[5], "\n\n"
                     , "p*=", p_opt_figS3[4], "\n\n"
                     , "p*=", p_opt_figS3[3], "\n\n"
                     , "p*=", p_opt_figS3[2], "\n\n"
                     , "p*=", p_opt_figS3[1], sep = ""
)
mtext(text = p_star_vals, 
      side = 4, line = 1, las=1, cex=0.8, outer=FALSE)
# ROW 3: B, varying over aw -------------------
N_figS3 <- array(NaN, c(no_states, resolution_B))
Y_figS3 <- array(NaN, c(no_states, resolution_B))
p_opt_figS3 <- array(NaN, c(no_states))
D_over_p <- D_B(p_vector, b)

for (i in 1:no_states) { # over parameter values
  aw <- aw_list[i]
  N_figS3[i,] <- N_eq14(p_vector, b, total_area, AAw, AAc, aw, ac)
  Y_figS3[i,] <- yield_B(p_vector, b, total_area, AAw, AAc, aw, ac)
  p_opt_figS3[i] <- p_vector[which(Y_figS3[i,]==max(Y_figS3[i,]))]
}
aw <- default_aw

plot(NULL
     , ylim= c(0, Ylim_B)
     , xlim= c(0, 1)
     , xaxt= 'n'
     #, yaxt= 'n'
)
for (i in 1:no_states) { # over parameter values
  lines(Y_figS3[i,]~D_B(p_vector, b)
        , col = F_colors[i]
        , lwd = linewidth
        , lty = F_lines
  )
  lines(N_figS3[i,]~D_B(p_vector, b)
        , col = N_colors[i]
        , lwd = linewidth
        , lty = N_lines
  )
  abline(v=D_B(p_opt_figS3[i], b), lty = opt_w_lines, lwd=2, col=F_colors[i])
}
# Letter label for this panel
text(offsetx, Ylim_B*offsety, "e", cex=2, col = "black")
# Legends & Row titles in left margin -----------------------------
text(leg_x, Ylim_B, "Varying aw \n(shape of N(w))"
     , cex=1, col="black", pos=4, xpd=NA)
legend(x=leg_x
       , y=Ylim_B*leg_y
       , xpd = NA
       , legend = paste("c=", c_states)
       , text.col = "transparent"
       , col = N_colors
       , pch = 19
       , cex = 0.9
       , pt.cex = 2
       , bty = 'n'
       , bg= "transparent"
)
legend(x=leg_x+0.1
       , y=Ylim_B*leg_y
       , xpd = NA
       , legend = paste("aw=", aw_list)
       , col = F_colors
       , pch = 19
       , cex = 0.9
       , pt.cex = 2
       , bty = 'n'
       , bg= "transparent"
)
# 2nd COLUMN ROW 3 -------------------------------
hatched_barplot_horizontal(rep(xbars, no_states))
barplot((1-p_opt_figS3)~b_vector
        , col = F_colors
        , xlim = c(0, xbars)
        , horiz = T
        , yaxt = 'n'
        , xaxt = 'n'
        , xpd = FALSE
        , add = TRUE
        , lwd = 2
)
box(lwd = 1)
# Letter label for this panel
text(xbars*(1-offsetx), 5.5, "f", cex=2, col = "black")
p_opt_figS3 <- round(p_opt_figS3, 2)
p_star_vals <- paste("p*=", p_opt_figS3[5], "\n\n"
                     , "p*=", p_opt_figS3[4], "\n\n"
                     , "p*=", p_opt_figS3[3], "\n\n"
                     , "p*=", p_opt_figS3[2], "\n\n"
                     , "p*=", p_opt_figS3[1], sep = ""
)
mtext(text = p_star_vals, 
      side = 4, line = 1, las=1, cex=0.8, outer=FALSE)

# ROW 4: B, varying over ac -------------------
N_figS3 <- array(NaN, c(no_states, resolution_B))
Y_figS3 <- array(NaN, c(no_states, resolution_B))
p_opt_figS3 <- array(NaN, c(no_states))
D_over_p <- D_B(p_vector, b)

for (i in 1:no_states) { # over parameter values
  ac <- ac_list[i]
  N_figS3[i,] <- N_eq14(p_vector, b, total_area, AAw, AAc, aw, ac)
  Y_figS3[i,] <- yield_B(p_vector, b, total_area, AAw, AAc, aw, ac)
  p_opt_figS3[i] <- p_vector[which(Y_figS3[i,]==max(Y_figS3[i,]))]
}
ac <- default_ac

plot(NULL
     , ylim= c(0, Ylim_B)
     , xlim= c(0, 1)
     #, yaxt= 'n'
)
for (i in 1:no_states) { # over parameter values
  lines(Y_figS3[i,]~D_B(p_vector, b)
        , col = F_colors[i]
        , lwd = linewidth
        , lty = F_lines
  )
  lines(N_figS3[i,]~D_B(p_vector, b)
        , col = N_colors[i]
        , lwd = linewidth
        , lty = N_lines
  )
  abline(v=D_B(p_opt_figS3[i], b), lty = opt_w_lines, lwd=2, col=F_colors[i])
}
mtext(text = "Proportion of pollinators\nin patch on crop: D = c/(c+w)", 
      side = 1, line = 4, las=1, cex=0.8, outer=FALSE)
# Letter label for this panel
text(offsetx, Ylim_B*offsety, "g", cex=2, col = "black")
# Legends in left margin -----------------
text(leg_x, Ylim_B, "Varying ac \n(shape of N(c))"
     , cex=1, col="black", pos=4, xpd=NA)
legend(x=leg_x
       , y=Ylim_B*leg_y
       , xpd = NA
       , legend = paste("c=", c_states)
       , text.col = "transparent"
       , col = N_colors
       , pch = 19
       , cex = 0.9
       , pt.cex = 2
       , bty = 'n'
       , bg= "transparent"
)
legend(x=leg_x+0.1
       , y=Ylim_B*leg_y
       , xpd = NA
       , legend = paste("ac=", ac_list)
       , col = F_colors
       , pch = 19
       , cex = 0.9
       , pt.cex = 2
       , bty = 'n'
       , bg= "transparent"
)
legend(x=leg_x
       , y=Ylim_B*0.2
       , xpd = NA
       , legend = c("N: Total bees \nin patch", "Y: 'yield'\n=ND(1-p)")
       , col = c(N_colors[3],F_colors[3])
       , lwd = 2
       , lty = c(N_lines, F_lines)
       , cex = 1
       , bty = 'n'
       , bg= "transparent"
       , y.intersp = 1.5
       )
# 2nd COLUMN ROW 4 -------------------------------
hatched_barplot_horizontal(rep(xbars, no_states))
barplot((1-p_opt_figS3)~b_vector
        , col = F_colors
        , xlim = c(0, xbars)
        , horiz = T
        , yaxt = 'n'
        , add = TRUE
)
box(lwd = 1)
mtext(text = "Fraction of total area that \nis crop to maximize Y (1-p*)", 
      side = 1, line = 4, las=1, cex=0.8, outer=FALSE)
# Letter label for this panel
text(xbars*(1-offsetx), 5.5, "h", cex=2, col = "black")
p_opt_figS3 <- round(p_opt_figS3, 2)
p_star_vals <- paste("p*=", p_opt_figS3[5], "\n\n"
                     , "p*=", p_opt_figS3[4], "\n\n"
                     , "p*=", p_opt_figS3[3], "\n\n"
                     , "p*=", p_opt_figS3[2], "\n\n"
                     , "p*=", p_opt_figS3[1], sep = ""
)
mtext(text = p_star_vals, 
      side = 4, line = 1, las=1, cex=0.8, outer=FALSE)

### END FIG S3 --------------------------



## ALTERNATIVE FIGURES


### FIG S3alt PARAMETERS ON p* 2x2 ------------------------
# x-axis D, y-axis y, 3 columns = 3 models
# 3 rows: Aw, aw, ac in colors, M and qa, qc for model A.3
# Figure format ---------------------
pointsize <- 0.6
offsetx <- 0.08
offsety <- 0.9
point_transparency <- 0.5
linewidth <- 4
linetransp <- 0.8
p_max <- 0.5
Ylim_B <- 1.5 # Max yield for graph
Y_clim_B <- 2 # Max for yield/c graph

# Make 4x3 panels
par(mfcol=c(2,2))
# Various other format adjustments (margins in order bottom, left, top, right)
par(oma = c(4,4,0,0), mar=c(0.5, 0.5, 0.5, 0.5), mgp=c(3, 1, 0), las=1)
# Parameters --------------------
default_b <- b
default_AAw <- AAw
default_aw <- aw
default_ac <- ac

# PANEL A: B, varying over b -------------------
N_figS3 <- array(NaN, c(no_states, resolution_B))
Y_figS3 <- array(NaN, c(no_states, resolution_B))
p_opt_figS3 <- array(NaN, c(no_states))

for (i in 1:no_states) { # over parameter values
  b <- b_vector[i]
  D_over_p <- D_B(p_vector, b)
  N_figS3[i,] <- N_eq14(p_vector, b, total_area, AAw, AAc, aw, ac)
  Y_figS3[i,] <- yield_B(p_vector, b, total_area, AAw, AAc, aw, ac)
  p_opt_figS3[i] <- p_vector[which( Y_figS3[i,]==max(Y_figS3[i,]))]
}
b <- default_b
D_over_p <- D_B(p_vector, b)

plot(NULL
     , ylim= c(0, Ylim_B)
     , xlim= c(0, 1)
     , xaxt= 'n'
)
for (i in 1:no_states) { # over parameter values
  lines(Y_figS3[i,]~D_B(p_vector, b_vector[i])
        , col = par_colors[i]
        , lwd = linewidth
  )
  abline(v=D_B(p_opt_figS3[i], b_vector[i]), lty = 2, lwd=2, col=par_colors_dark[i])
}
# Letter label for this panel
text(0.5, Ylim_B*offsety, "a", cex=2, col = "black")
legend("topleft"
       , legend = paste("b=", b_vector)
       , col = par_colors_dark
       , pch = 19
       , cex = 0.9
       , pt.cex = 2
       , bty = 'n'
       , bg= "transparent"
)
# PANEL B: B, varying over AAw -------------------
N_figS3 <- array(NaN, c(no_states, resolution_B))
Y_figS3 <- array(NaN, c(no_states, resolution_B))
p_opt_figS3 <- array(NaN, c(no_states))
D_over_p <- D_B(p_vector, b)

for (i in 1:no_states) { # over parameter values
  AAw <- AAw_list[i]
  N_figS3[i,] <- N_eq14(p_vector, b, total_area, AAw, AAc, aw, ac)
  Y_figS3[i,] <- yield_B(p_vector, b, total_area, AAw, AAc, aw, ac)
  p_opt_figS3[i] <- p_vector[which(Y_figS3[i,]==max(Y_figS3[i,]))]
}
AAw <- default_AAw

plot(NULL
     , ylim= c(0, Ylim_B)
     , xlim= c(0, 1)
)
for (i in 1:no_states) { # over parameter values
  lines(Y_figS3[i,]~D_B(p_vector, b)
        , col = par_colors[i]
        , lwd = linewidth
  )
  abline(v=D_B(p_opt_figS3[i], b), lty = 2, lwd=2, col=par_colors_dark[i])
}
# Letter label for this panel
text(offsetx, Ylim_B*offsety, "b", cex=2, col = "black")
legend("top"
       , legend = paste("A=", AAw_list)
       , col = par_colors_dark
       , pch = 19
       , cex = 0.9
       , pt.cex = 2
       , bty = 'n'
       , bg= "transparent"
)
# PANEL C: B, varying over aw -------------------
N_figS3 <- array(NaN, c(no_states, resolution_B))
Y_figS3 <- array(NaN, c(no_states, resolution_B))
p_opt_figS3 <- array(NaN, c(no_states))
D_over_p <- D_B(p_vector, b)

for (i in 1:no_states) { # over parameter values
  aw <- aw_list[i]
  N_figS3[i,] <- N_eq14(p_vector, b, total_area, AAw, AAc, aw, ac)
  Y_figS3[i,] <- yield_B(p_vector, b, total_area, AAw, AAc, aw, ac)
  p_opt_figS3[i] <- p_vector[which(Y_figS3[i,]==max(Y_figS3[i,]))]
}
aw <- default_aw

plot(NULL
     , ylim= c(0, Ylim_B)
     , xlim= c(0, 1)
     , xaxt= 'n'
     , yaxt= 'n'
)
for (i in 1:no_states) { # over parameter values
  lines(Y_figS3[i,]~D_B(p_vector, b)
        , col = par_colors[i]
        , lwd = linewidth
  )
  abline(v=D_B(p_opt_figS3[i], b), lty = 2, lwd=2, col=par_colors_dark[i])
}
# Letter label for this panel
text(0.5, Ylim_B*offsety, "c", cex=2, col = "black")
legend("topleft"
       , legend = paste("aw=", aw_list)
       , col = par_colors_dark
       , pch = 19
       , cex = 0.9
       , pt.cex = 2
       , bty = 'n'
       , bg= "transparent"
)
# PANEL D: B, varying over ac -------------------
N_figS3 <- array(NaN, c(no_states, resolution_B))
Y_figS3 <- array(NaN, c(no_states, resolution_B))
p_opt_figS3 <- array(NaN, c(no_states))
D_over_p <- D_B(p_vector, b)

for (i in 1:no_states) { # over parameter values
  ac <- ac_list[i]
  N_figS3[i,] <- N_eq14(p_vector, b, total_area, AAw, AAc, aw, ac)
  Y_figS3[i,] <- yield_B(p_vector, b, total_area, AAw, AAc, aw, ac)
  p_opt_figS3[i] <- p_vector[which(Y_figS3[i,]==max(Y_figS3[i,]))]
}
ac <- default_ac

plot(NULL
     , ylim= c(0, Ylim_B)
     , xlim= c(0, 1)
     , yaxt= 'n'
)
for (i in 1:no_states) { # over parameter values
  lines(Y_figS3[i,]~D_B(p_vector, b)
        , col = par_colors[i]
        , lwd = linewidth
  )
  abline(v=D_B(p_opt_figS3[i], b), lty = 2, lwd=2, col=par_colors_dark[i])
}
# Letter label for this panel
text(0.5, Ylim_B*offsety, "d", cex=2, col = "black")
legend("topleft"
       , legend = paste("ac=", ac_list)
       , col = par_colors_dark
       , pch = 19
       , cex = 0.9
       , pt.cex = 2
       , bty = 'n'
       , bg= "transparent"
)
# X-Axis label ----------------------------
# Text in outer margin for overarching x-axis label
mtext(text = "Proportion of pollinators in patch on crop: D = c/(c+w)", 
      side = 1, line = 2, las=1, cex=0.8, outer=TRUE)
# Y-Axis label ----------------------------
# Text in outer margin for overarching x-axis label
mtext(text = "'Yield' Y = Area of crop (p) * Bees on crop (N)", 
      side = 2, line = 2, las=3, cex=0.8, outer=TRUE)

### END FIG S3alt --------------------------

### FIG S3alt2 PARAMETERS ON p* 2x2 ------------------------
# x-axis D, y-axis y, 3 columns = 3 models
# 3 rows: Aw, aw, ac in colors, M and qa, qc for model A.3
# Figure format ---------------------
pointsize <- 0.6
offsetx <- 0.12
offsety <- 0.9
point_transparency <- 0.5
linewidth <- 4
linetransp <- 0.8
p_max <- 0.5
Ylim_B <- 2 # Max yield for graph
Y_clim_B <- 2 # Max for yield/c graph
leg_x <- -0.9

# Make 4x3 panels
par(mfcol=c(4,1))
# Various other format adjustments (margins in order bottom, left, top, right)
par(oma = c(5,12,0.5,0), mar=c(0.5, 0.5, 0.5, 0.5), mgp=c(3, 1, 0), las=1)
# Parameters --------------------
default_b <- b
default_AAw <- AAw
default_aw <- aw
default_ac <- ac

# PANEL A: B, varying over b -------------------
N_figS3 <- array(NaN, c(no_states, resolution_B))
Y_figS3 <- array(NaN, c(no_states, resolution_B))
p_opt_figS3 <- array(NaN, c(no_states))

for (i in 1:no_states) { # over parameter values
  b <- b_vector[i]
  D_over_p <- D_B(p_vector, b)
  N_figS3[i,] <- N_eq14(p_vector, b, total_area, AAw, AAc, aw, ac)
  Y_figS3[i,] <- yield_B(p_vector, b, total_area, AAw, AAc, aw, ac)
  p_opt_figS3[i] <- p_vector[which( Y_figS3[i,]==max(Y_figS3[i,]))]
}
b <- default_b
D_over_p <- D_B(p_vector, b)

plot(NULL
     , ylim= c(0, Ylim_B)
     , xlim= c(0, 1)
     , xaxt= 'n'
)
for (i in 1:no_states) { # over parameter values
  lines(Y_figS3[i,]~D_B(p_vector, b_vector[i])
        , col = F_colors[i]
        , lwd = linewidth
        , lty = F_lines
  )
  lines(N_figS3[i,]~D_B(p_vector, b_vector[i])
        , col = N_colors[i]
        , lwd = linewidth
        , lty = N_lines
  )
  abline(v=D_B(p_opt_figS3[i], b_vector[i]), lty = opt_w_lines, lwd=2, col=F_colors[i])
}
# Letter label for this panel
text(offsetx, Ylim_B*offsety, "a", cex=2, col = "black")
# Legends & Row titles in left margin -----------------------------
text(leg_x, Ylim_B, "Varying c \n(amount of crop)"
     , cex=1, col="black", pos=4, xpd=NA)
legend(x=leg_x
       , y=Ylim_B*leg_y
       , xpd = NA
       , legend = paste("c=", c_states)
       , text.col = "transparent"
       , col = N_colors
       , pch = 19
       , cex = 0.9
       , pt.cex = 2
       , bty = 'n'
       , bg= "transparent"
)
legend(x=leg_x+0.1
       , y=Ylim_B*leg_y
       , xpd = NA
       , legend = paste("b=", b_vector)
       , col = F_colors
       , pch = 19
       , cex = 0.9
       , pt.cex = 2
       , bty = 'n'
       , bg= "transparent"
)
# PANEL B: B, varying over AAw -------------------
N_figS3 <- array(NaN, c(no_states, resolution_B))
Y_figS3 <- array(NaN, c(no_states, resolution_B))
p_opt_figS3 <- array(NaN, c(no_states))
D_over_p <- D_B(p_vector, b)

for (i in 1:no_states) { # over parameter values
  AAw <- AAw_list[i]
  N_figS3[i,] <- N_eq14(p_vector, b, total_area, AAw, AAc, aw, ac)
  Y_figS3[i,] <- yield_B(p_vector, b, total_area, AAw, AAc, aw, ac)
  p_opt_figS3[i] <- p_vector[which(Y_figS3[i,]==max(Y_figS3[i,]))]
}
AAw <- default_AAw

plot(NULL
     , ylim= c(0, Ylim_B)
     , xlim= c(0, 1)
     , xaxt='n'
)
for (i in 1:no_states) { # over parameter values
  lines(Y_figS3[i,]~D_B(p_vector, b)
        , col = F_colors[i]
        , lwd = linewidth
        , lty = F_lines
  )
  lines(N_figS3[i,]~D_B(p_vector, b)
        , col = N_colors[i]
        , lwd = linewidth
        , lty = N_lines
  )
  abline(v=D_B(p_opt_figS3[i], b), lty = opt_w_lines, lwd=2, col=F_colors[i])
}
# Letter label for this panel
text(offsetx, Ylim_B*offsety, "b", cex=2, col = "black")
# Legends & Row titles in left margin -----------------------------
text(leg_x, Ylim_B, "Varying Aw \n(attractivenes of w)"
     , cex=1, col="black", pos=4, xpd=NA)
legend(x=leg_x
       , y=Ylim_B*leg_y
       , xpd = NA
       , legend = paste("c=", c_states)
       , text.col = "transparent"
       , col = N_colors
       , pch = 19
       , cex = 0.9
       , pt.cex = 2
       , bty = 'n'
       , bg= "transparent"
)
legend(x=leg_x+0.1
       , y=Ylim_B*leg_y
       , xpd = NA
       , legend = paste("Aw=", AAw_list)
       , col = F_colors
       , pch = 19
       , cex = 0.9
       , pt.cex = 2
       , bty = 'n'
       , bg= "transparent"
)
# PANEL C: B, varying over aw -------------------
N_figS3 <- array(NaN, c(no_states, resolution_B))
Y_figS3 <- array(NaN, c(no_states, resolution_B))
p_opt_figS3 <- array(NaN, c(no_states))
D_over_p <- D_B(p_vector, b)

for (i in 1:no_states) { # over parameter values
  aw <- aw_list[i]
  N_figS3[i,] <- N_eq14(p_vector, b, total_area, AAw, AAc, aw, ac)
  Y_figS3[i,] <- yield_B(p_vector, b, total_area, AAw, AAc, aw, ac)
  p_opt_figS3[i] <- p_vector[which(Y_figS3[i,]==max(Y_figS3[i,]))]
}
aw <- default_aw

plot(NULL
     , ylim= c(0, Ylim_B)
     , xlim= c(0, 1)
     , xaxt= 'n'
     #, yaxt= 'n'
)
for (i in 1:no_states) { # over parameter values
  lines(Y_figS3[i,]~D_B(p_vector, b)
        , col = F_colors[i]
        , lwd = linewidth
        , lty = F_lines
  )
  lines(N_figS3[i,]~D_B(p_vector, b)
        , col = N_colors[i]
        , lwd = linewidth
        , lty = N_lines
  )
  abline(v=D_B(p_opt_figS3[i], b), lty = opt_w_lines, lwd=2, col=F_colors[i])
}
# Letter label for this panel
text(offsetx, Ylim_B*offsety, "c", cex=2, col = "black")
# Legends & Row titles in left margin -----------------------------
text(leg_x, Ylim_B, "Varying aw \n(shape of N(w))"
     , cex=1, col="black", pos=4, xpd=NA)
legend(x=leg_x
       , y=Ylim_B*leg_y
       , xpd = NA
       , legend = paste("c=", c_states)
       , text.col = "transparent"
       , col = N_colors
       , pch = 19
       , cex = 0.9
       , pt.cex = 2
       , bty = 'n'
       , bg= "transparent"
)
legend(x=leg_x+0.1
       , y=Ylim_B*leg_y
       , xpd = NA
       , legend = paste("aw=", aw_list)
       , col = F_colors
       , pch = 19
       , cex = 0.9
       , pt.cex = 2
       , bty = 'n'
       , bg= "transparent"
)
# PANEL D: B, varying over ac -------------------
N_figS3 <- array(NaN, c(no_states, resolution_B))
Y_figS3 <- array(NaN, c(no_states, resolution_B))
p_opt_figS3 <- array(NaN, c(no_states))
D_over_p <- D_B(p_vector, b)

for (i in 1:no_states) { # over parameter values
  ac <- ac_list[i]
  N_figS3[i,] <- N_eq14(p_vector, b, total_area, AAw, AAc, aw, ac)
  Y_figS3[i,] <- yield_B(p_vector, b, total_area, AAw, AAc, aw, ac)
  p_opt_figS3[i] <- p_vector[which(Y_figS3[i,]==max(Y_figS3[i,]))]
}
ac <- default_ac

plot(NULL
     , ylim= c(0, Ylim_B)
     , xlim= c(0, 1)
     #, yaxt= 'n'
)
for (i in 1:no_states) { # over parameter values
  lines(Y_figS3[i,]~D_B(p_vector, b)
        , col = F_colors[i]
        , lwd = linewidth
        , lty = F_lines
  )
  lines(N_figS3[i,]~D_B(p_vector, b)
        , col = N_colors[i]
        , lwd = linewidth
        , lty = N_lines
  )
  abline(v=D_B(p_opt_figS3[i], b), lty = opt_w_lines, lwd=2, col=F_colors[i])
}
# Letter label for this panel
text(offsetx, Ylim_B*offsety, "d", cex=2, col = "black")
# Legends in left margin -----------------
text(leg_x, Ylim_B, "Varying ac \n(shape of N(c))"
     , cex=1, col="black", pos=4, xpd=NA)
legend(x=leg_x
       , y=Ylim_B*leg_y
       , xpd = NA
       , legend = paste("c=", c_states)
       , text.col = "transparent"
       , col = N_colors
       , pch = 19
       , cex = 0.9
       , pt.cex = 2
       , bty = 'n'
       , bg= "transparent"
)
legend(x=leg_x+0.1
       , y=Ylim_B*leg_y
       , xpd = NA
       , legend = paste("ac=", ac_list)
       , col = F_colors
       , pch = 19
       , cex = 0.9
       , pt.cex = 2
       , bty = 'n'
       , bg= "transparent"
)
legend(x=-1.
       , y=Ylim_B*0.2
       , xpd = NA
       , legend = c("N: Total bees \nin patch", "Y: 'yield'\n=ND(1-p)")
       , col = c(N_colors[3],F_colors[3])
       , lwd = 2
       , lty = c(N_lines, F_lines)
       , cex = 1
       , bty = 'n'
       , bg= "transparent"
       , y.intersp = 1.5
)
# X-Axis label ----------------------------
# Text in outer margin for overarching x-axis label
mtext(text = "Proportion of pollinators\nin patch on crop: D = c/(c+w)", 
      side = 1, line = 3, las=1, cex=0.8, outer=TRUE)
# Y-Axis label ----------------------------
# Text in outer margin for overarching x-axis label
#mtext(text = "'Yield' Y = Area of crop (p) * Bees on crop (N)", 
#      side = 2, line = 2, las=3, cex=0.8, outer=TRUE)
### END FIG S3alt2 --------------------------

