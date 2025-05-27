# Alasdair Houston & Anna Dornhaus
# R script for numerical model analysis and figures
# Spillover project: Bees 'spilling over' into crop from wildflowers planted to attract them

# Paper reference: XXX (to be updated on acceptance of the paper)
### PAPER VERSION R1 April 2025
# FIGURES IN PAPER

### Libraries and graphics setup ------------------
# Libraries purely for graphics
library(scales)
library(viridis)
### COLORS 1
# colorlist has the colors for N, Fc, Nc, D, Fw, Nw, the tangent in this order
templist1 <- inferno(6)
templist2 <- mako(6)
colorlist <- c(templist1[2], templist1[3], templist1[4], templist1[5], templist2[4], templist2[5], templist2[3])

### COLORS 2
# colors for parameter variations
no_states <- 7
par_colors_dark <- viridis(no_states)
par_colors <- alpha(par_colors_dark, 0.5)

# Resetting graphics parameters
#dev.off()
opar <- par()


### GENERAL PARAMETERS ------------------

# Default parameter set
# Model A.1 
aw <- 0.4 # diminishing returns of attraction based on 'amount' of flowers 
AAw <- 10   # constant translating 'amount' of wildflowers into pollinators
AAc <- 1 # constant translating 'amount' of crop into pollinators
ac <- 0.4
# Fig. 2
l <- 10
rho <- 0.02
gamma <- 3
 
## For numerical calculations
max_w <- 1000 # for numerical calculations
resolution <- 30000 # number of calculated values
w_list <- seq(from=0, to=max_w, length.out=resolution)
c_list <- w_list

## For plotting
pointsize <- 0.6
Flim <- 500 # Max F for graph
F_clim <- Flim / 10 # Max for F/c graph
wlim <- 50 # for plotting
clim <- 250 # Max c for graph
offsetx <- 0.80
offsety <- 0.98


###### MODEL A ------------------
# Functions -----------------------

# General for all models
# equation for D
D <- function(c,w) {
  return(c/(c+w))
}
# equation for F
F <- function(N, D) {
  return(N*D)
}

# Equation [5] defines a possible N(w); not assuming anything about Nc
N_eq5 <- function(w, Nc, a, A) {
  return(Nc + (A*w) / (a+w))
}

# This is our equation for optimal w as a function of c and Nc
wopt_eq8 <- function(c, Nc, a, A) {
  # Naming the parameters of the quadratic qa, qb, and qc here
  qa <- A + Nc
  qb <- 2*a*Nc
  qc <- Nc*a^2 - A*a*c
  w_opt <- (-qb + sqrt(qb^2 - 4*qa*qc)) / (2*qa)
  return(w_opt)
}

# General function of N(c,w)
Nof_candw <- function(c,w) {
  return((AAw*w) / (aw+w) + (AAc*c) / (ac+c))
}

# Function for drawing the tangent in Fig. 1
tangent_N <- function(x, c_fixed, Nc) {
  slope <- Nof_candw(c_fixed,wopt_eq8(c_fixed, Nc,aw, AAw))/(c_fixed+wopt_eq8(c_fixed, Nc,aw, AAw))
  return((x+c_fixed)*slope)  
}


### ACTUAL MANUSCRIPT FIGURES ------------------------------

# Fig. 1 Conceptual argument & graphical solution ----------
# All par() functions set margins and other details about figures.
par(mfrow=c(1,1))
# Margins are set in order bottom, left, top, right
par(mar=c(0.5, 0.5, 0.5, 0.5), oma=c(4, 0, 0, 0), mgp=c(3, 1, 0), las=1)

# PANEL A
# Graphical solution with tangent
# We show, for a fixed c, how N(w), F(w), and D(w) behave.
c_fixed <- 8
Nc <- Nof_candw(c_fixed, 0)
w_optimal <- wopt_eq8(c_fixed, Nc, aw, AAw)
ymax <- 25

# First set up the frame of the plot  
plot(NULL
     , xlab = "" 
     , ylab ="Number of bees"
     , ylim = c(0, ymax)
     , xlim = c(-10, 15)
     , yaxt = 'n'
)
# X-axis label with odd spacing
mtext(text = "            -c (amount of crop)                             w (amount of wildflowers)",
      side = 1, line = 3, adj = 0, cex = 1)
# N(w) according to eq [5]; 'x' is w
curve(N_eq5(x, Nc, aw, AAw)
      , from=0
      , to=15
      , add=TRUE
      , lty=1
      , lwd=3
      , col=colorlist[1]
      , n = 1000
)
# Tangent line dashed
curve(tangent_N(x, c_fixed, Nc)
      , add=TRUE
      , lty=3
      , lwd=2
      , col=colorlist[7]
      )
# Or as line segment just from x-axis to tangent point
lines(c(-c_fixed, w_optimal), c(0, N_eq5(w_optimal, Nc, aw, AAw))
      , lty=1
      , lwd=3
      , col=colorlist[7]
      )
# F
curve(F(N_eq5(x, Nc, aw, AAw), D(c_fixed, x))
      , from = 0
      , to = 15
      , add=TRUE
      , lty=1
      , lwd=3
      , col=colorlist[2]
      , n = 1000
      )
abline(h=0, lty=3) # Horizontal line through origin
abline(v=0, lty=3) # Vertical line through origin
# Drawing a horizontal line segment at Nc
segments(x0 = 0, x1 = 15, y0 = Nc, y1 = Nc, lty=1, lwd=3, col=colorlist[3]) 
# thin vertical lines at c and w*
abline(v=w_optimal, lty=2, lwd=2, col="grey")
abline(v=-c_fixed, lty=2, lwd=2, col="grey")
# Labels for c and w*
text(1.5-c_fixed, ymax*0.6, paste("c=", round(c_fixed, 1)), cex=1.5)
text(1.5+w_optimal, ymax*0.6, paste("w*=", round(w_optimal,2)), cex=1.5)
# Overall legend
legend("topleft"
       , legend = c("N (bees on c+w)", "F (bees on crop with w)", "Nc (bees on crop without w)", "Tangent from -c")
       , col = c(colorlist[1], colorlist[2], colorlist[3], colorlist[7])
       , pch = 19
       , cex = 1.2
       , bg="transparent"
)


# Fig. 2 Different behaviors of N --------------------------------

# Setting parameters for the figure
c_fixed <- 8
Nc <- Nof_candw(c_fixed, 0)
w_optimal <- wopt_eq8(c_fixed, Nc, aw, AAw)
ymax <- 25
# Since we use a discrete list of c values, this determines which index is 
# closest to the actual c_fixed defined above
diff_to_cfixed <- abs(c_list-c_fixed)
c_fixed_index <- which(diff_to_cfixed == min(diff_to_cfixed))

# Basic, over w_list
D_over_w <- D(c_fixed, w_list)
# varying w, keeping c fixed, and calculating the resulting D
# (nothing optimized here, no w*)
# Same:
N_over_w <- N_eq5(w_list, Nc, aw, AAw)
F_over_w <- F(N_over_w, D_over_w)

# Same over a range of cs
D_over_c <- D(c_list, w_optimal)
N_over_c <- N_eq5(w_optimal, AAc*(c_list), aw, AAw)
F_over_c <- F(N_over_c, D_over_c)

# par() sets margins (in order bottom, left, top, right) and other graphical
# parameters
par(mfrow=c(3,2)) # 3 rows and 2 columns of graphs
par(mar=c(0.5, 0.5, 0.5, 0.5), oma=c(4, 5, 0, 0), mgp=c(3, 1, 0), las=1)
# Setting parameters to place text consistently
textposx <- 50
textposy <- 0.35

# First part:
# Quick numerical illustration of N vs D
# N1: N equal to c+w
plot(D_over_w ~ w_list
     , col = colorlist[4]
     , pch = 1
     , xlim = c(0, wlim)
     , ylim = c(0,1)
     , xlab = "w (amount of wildflowers)"
     , ylab = "Effect on crop"
     , xaxt = 'n'
     , cex = pointsize
     , yaxp = c(0, 2, 2)
)
N1 <- function(c, w) {
  return(AAw*(c+w))
}
# N/N(wlim)
stand_N <- N1(c_fixed, w_list)/N1(c_fixed, wlim) 
# Fc/N(wlim)
Fc <- F(D_over_w, N1(c_fixed, w_list))
stand_Fc <- Fc/N1(c_fixed, wlim)
polygon(c(w_list, rev(w_list)), c(c(1:resolution)*0, rev(stand_Fc))
        , col = colorlist[2]
        , lty = 0
)
# Fw
polygon(c(w_list, rev(w_list)), c(stand_Fc, rev(stand_N))
        , col = colorlist[5]
        , lty = 0
)
#Fw <- (N1(c_fixed, w_list) - F(D_over_w, N1(c_fixed, w_list)))
#stand_Fw <- Fw/N1(c_fixed, wlim)
# Fc
points(stand_Fc ~ w_list
       , col = colorlist[2]
       , pch = 19
       , cex = pointsize
)
# Nw
Nw <- N1(0, w_list)/N1(c_fixed, wlim)
points(Nw ~ w_list
       , col = colorlist[6]
       , pch = 1
       , cex = pointsize
)
# N/N(wlim)
points(stand_N ~ w_list
       , col = colorlist[1]
       , pch = 1
       , cex = pointsize
)
# Nc
Nc <- N1(c_fixed, 0)/N1(c_fixed, wlim)
abline(h=Nc, lty=2, lwd=2, col=colorlist[3])
# D
points(D_over_w ~ w_list
       , col = colorlist[4]
       , pch = 1
       , cex = pointsize
)
text(textposx, textposy, "N additive", cex=1.5, adj=c(0,0), pos=2)
text(7, 0.9, "a", cex=2, adj=c(0,0))
legend("topright"
       , legend = c("N (total bees)", "F (bees on c)", "Nc (bees on c without w)", "D (fraction on c)", "(bees on w)", "(bees on w without c)")
       , col = c(colorlist[1], colorlist[2], colorlist[3], colorlist[4], colorlist[5], colorlist[6])
       , pch = 19
       , cex = 1.2
       #       , box.lty = 0
)

# N2: N equal to c + w + l
plot(D_over_w ~ w_list
     , col = colorlist[4]
     , pch = 19
     , xlim = c(0, wlim)
     , ylim = c(0,1)
     , xlab = "w (amount of wildflowers)"
     , ylab = "Effect on crop"
     , xaxt = 'n'
     , yaxt = 'n'
     , cex = pointsize
     , yaxp = c(0, 2, 2)
)
N2 <- function(c, w) {
  return(w + c + l)
}
# N/N(wlim)
stand_N <- N2(c_fixed, w_list)/N2(c_fixed, wlim) 
# Fc/N(wlim)
Fc <- F(D_over_w, N2(c_fixed, w_list))
stand_Fc <- Fc/N2(c_fixed, wlim)
polygon(c(w_list, rev(w_list)), c(c(1:resolution)*0, rev(stand_Fc))
        #        , data = cw_ModelB
        , col = colorlist[2]
        , lty = 0
)
# Fw
#Fw <- (N1(c_fixed, w_list) - F(D_over_w, N1(c_fixed, w_list)))
#stand_Fw <- Fw/N1(c_fixed, wlim)
polygon(c(w_list, rev(w_list)), c(stand_Fc, rev(stand_N))
        , col = colorlist[5]
        , lty = 0
)
# Fc
points(stand_Fc ~ w_list
       , col = colorlist[2]
       , pch = 19
       , cex = pointsize
)
# Nw
Nw <- N2(0, w_list)/N2(c_fixed, wlim)
points(Nw ~ w_list
       , col = colorlist[6]
       , pch = 1
       , cex = pointsize
)
# N/N(wlim)
points(stand_N ~ w_list
       , col = colorlist[1]
       , pch = 1
       , cex = pointsize
)
# Nc
Nc <- N2(c_fixed, 0)/N2(c_fixed, wlim)
abline(h=Nc, lty=2, lwd=2, col=colorlist[3])
# D
points(D_over_w ~ w_list
       , col = colorlist[4]
       , pch = 1
       , cex = pointsize
)
text(textposx, textposy, "N with local bees", cex=1.5, adj=c(0,0), pos=2)
text(7, 0.9, "b", cex=2, adj=c(0,0))

# N3: N equal to c + w * A
plot(D_over_w ~ w_list
     , col = colorlist[4]
     , pch = 19
     , xlim = c(0, wlim)
     , ylim = c(0,1)
     , xlab = "w (amount of wildflowers)"
     , ylab = "Effect on crop"
     , xaxt = 'n'
     #, yaxt = 'n'
     , cex = pointsize
     , yaxp = c(0, 2, 2)
)
N3 <- function(c, w) {
  return(AAw*w + c)
}
# N/N(wlim)
stand_N <- N3(c_fixed, w_list)/N3(c_fixed, wlim) 
# Fc/N(wlim)
Fc <- F(D_over_w, N3(c_fixed, w_list))
stand_Fc <- Fc/N3(c_fixed, wlim)
polygon(c(w_list, rev(w_list)), c(c(1:resolution)*0, rev(stand_Fc))
        #        , data = cw_ModelB
        , col = colorlist[2]
        , lty = 0
)
# Fw
#Fw <- (N1(c_fixed, w_list) - F(D_over_w, N1(c_fixed, w_list)))
#stand_Fw <- Fw/N1(c_fixed, wlim)
polygon(c(w_list, rev(w_list)), c(stand_Fc, rev(stand_N))
        , col = colorlist[5]
        , lty = 0
)
# Fc
points(stand_Fc ~ w_list
       , col = colorlist[2]
       , pch = 19
       , cex = pointsize
)
# Nw
Nw <- N3(0, w_list)/N3(c_fixed, wlim)
points(Nw ~ w_list
       , col = colorlist[6]
       , pch = 1
       , cex = pointsize
)
# N/N(wlim)
points(stand_N ~ w_list
       , col = colorlist[1]
       , pch = 1
       , cex = pointsize
)
# Nc
Nc <- N3(c_fixed, 0)/N3(c_fixed, wlim)
abline(h=Nc, lty=2, lwd=2, col=colorlist[3])
# D
points(D_over_w ~ w_list
       , col = colorlist[4]
       , pch = 1
       , cex = pointsize
)
text(textposx, textposy, "attractive w", cex=1.5, adj=c(0,0), pos=2)
text(7, 0.9, "c", cex=2, adj=c(0,0))

# N4: N diminishing with c+w
plot(D_over_w ~ w_list
     , col = colorlist[4]
     , pch = 19
     , xlim = c(0, wlim)
     , ylim = c(0,1)
     , xlab = "w (amount of wildflowers)"
     , ylab = "Effect on crop"
     , xaxt = 'n'
     , yaxt = 'n'
     , cex = pointsize
     , yaxp = c(0, 2, 2)
)
N4 <- function(c, w) {
  return(c^aw+w^aw)
}
# N/N(wlim)
stand_N <- N4(c_fixed, w_list)/N4(c_fixed, wlim) 
# Fc/N(wlim)
Fc <- F(D_over_w, N4(c_fixed, w_list))
stand_Fc <- Fc/N4(c_fixed, wlim)
polygon(c(w_list, rev(w_list)), c(c(1:resolution)*0, rev(stand_Fc))
        #        , data = cw_ModelB
        , col = colorlist[2]
        , lty = 0
)
# Fw
#Fw <- (N1(c_fixed, w_list) - F(D_over_w, N1(c_fixed, w_list)))
#stand_Fw <- Fw/N1(c_fixed, wlim)
polygon(c(w_list, rev(w_list)), c(stand_Fc, rev(stand_N))
        , col = colorlist[5]
        , lty = 0
)
# Fc
points(stand_Fc ~ w_list
       , col = colorlist[2]
       , pch = 19
       , cex = pointsize
)
# Nw
Nw <- N4(0, w_list)/N4(c_fixed, wlim)
points(Nw ~ w_list
       , col = colorlist[6]
       , pch = 1
       , cex = pointsize
)
# N/N(wlim)
points(stand_N ~ w_list
       , col = colorlist[1]
       , pch = 1
       , cex = pointsize
)
# Nc
Nc <- N4(c_fixed, 0)/N4(c_fixed, wlim)
abline(h=Nc, lty=2, lwd=2, col=colorlist[3])
# D
points(D_over_w ~ w_list
       , col = colorlist[4]
       , pch = 1
       , cex = pointsize
)
text(textposx, textposy, "N decelerating", cex=1.5, adj=c(0,0), pos=2)
text(7, 0.9, "d", cex=2, adj=c(0,0))

# N5: N saturating with c+w
# bottom, left, top, right
plot(D_over_w ~ w_list
     , col = colorlist[4]
     , pch = 19
     , xlim = c(0, wlim)
     , ylim = c(0,1)
     , xlab = "w (amount of wildflowers)"
     , ylab = "Effect on crop"
     , cex = pointsize
     , yaxp = c(0, 2, 2)
)
N5 <- function(c, w) {
  return((1 - exp(-aw*w)) + (1 - exp(-aw*c)))
}
# N/N(wlim)
stand_N <- N5(c_fixed, w_list)/N5(c_fixed, wlim) 
# Fc/N(wlim)
Fc <- F(D_over_w, N5(c_fixed, w_list))
stand_Fc <- Fc/N5(c_fixed, wlim)
polygon(c(w_list, rev(w_list)), c(c(1:resolution)*0, rev(stand_Fc))
        #        , data = cw_ModelB
        , col = colorlist[2]
        , lty = 0
)
# Fw
#Fw <- (N1(c_fixed, w_list) - F(D_over_w, N1(c_fixed, w_list)))
#stand_Fw <- Fw/N1(c_fixed, wlim)
polygon(c(w_list, rev(w_list)), c(stand_Fc, rev(stand_N))
        , col = colorlist[5]
        , lty = 0
)
# Fc
points(stand_Fc ~ w_list
       , col = colorlist[2]
       , pch = 19
       , cex = pointsize
)
# Nw
Nw <- N5(0, w_list)/N5(c_fixed, wlim)
points(Nw ~ w_list
       , col = colorlist[6]
       , pch = 1
       , cex = pointsize
)
# N/N(wlim)
points(stand_N ~ w_list
       , col = colorlist[1]
       , pch = 1
       , cex = pointsize
)
# Nc
Nc <- N5(c_fixed, 0)/N5(c_fixed, wlim)
abline(h=Nc, lty=2, lwd=2, col=colorlist[3])
# D
points(D_over_w ~ w_list
       , col = colorlist[4]
       , pch = 1
       , cex = pointsize
)
text(textposx, textposy, "N saturating", cex=1.5, adj=c(0,0), pos=2)
text(7, 0.9, "e", cex=2, adj=c(0,0))
mtext(text = "w (amount of wildflowers)",
      side = 1, line = 3, cex = 0.8)

# N6: Holling III, S curve initially accelerating
# bottom, left, top, right
plot(D_over_w ~ w_list
     , col = colorlist[4]
     , pch = 19
     , xlim = c(0, wlim)
     , ylim = c(0,1)
     , xlab = "w (amount of wildflowers)"
     , ylab = "Effect on crop"
     , cex = pointsize
     , yaxp = c(0, 2, 2)
)
N6 <- function(c, w) {
  return( ((rho*(c+w))^gamma)/(1+((rho*(c+w))^gamma)))
}
# N/N(wlim)
stand_N <- N6(c_fixed, w_list)/N6(c_fixed, wlim) 
# Fc/N(wlim)
Fc <- F(D_over_w, N6(c_fixed, w_list))
stand_Fc <- Fc/N6(c_fixed, wlim)
polygon(c(w_list, rev(w_list)), c(c(1:resolution)*0, rev(stand_Fc))
        #        , data = cw_ModelB
        , col = colorlist[2]
        , lty = 0
)
# Fw
#Fw <- (N1(c_fixed, w_list) - F(D_over_w, N1(c_fixed, w_list)))
#stand_Fw <- Fw/N1(c_fixed, wlim)
polygon(c(w_list, rev(w_list)), c(stand_Fc, rev(stand_N))
        , col = colorlist[5]
        , lty = 0
)
# Fc
points(stand_Fc ~ w_list
       , col = colorlist[2]
       , pch = 19
       , cex = pointsize
)
# Nw
Nw <- N6(0, w_list)/N6(c_fixed, wlim)
points(Nw ~ w_list
       , col = colorlist[6]
       , pch = 1
       , cex = pointsize
)
# N/N(wlim)
points(stand_N ~ w_list
       , col = colorlist[1]
       , pch = 1
       , cex = pointsize
)
# Nc
Nc <- N6(c_fixed, 0)/N6(c_fixed, wlim)
abline(h=Nc, lty=2, lwd=2, col=colorlist[3])
# D
points(D_over_w ~ w_list
       , col = colorlist[4]
       , pch = 1
       , cex = pointsize
)
text(textposx, textposy, "N Holling III", cex=1.5, adj=c(0,0), pos=2)
text(7, 0.9, "f", cex=2, adj=c(0,0))
mtext(text = "w (amount of wildflowers)",
      side = 1, line = 3, cex = 0.8)

mtext(text = "Distribution of pollinator visits (scaled)",
      side = 2, line = 2, cex = 1, las = 3, outer = TRUE)

###### MODEL B ---------------------------------------

# Functions ---------------------

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

# Parameters overall ----------------------
b <- 1
p_vector <- seq(0, 1, length.out=resolution)
b_vector <- seq(1,100, length.out=no_states)
total_area <- 1

# Parameters panel a ---------------------
AAw <- 1
AAc <- 1
Flim_B <- 1.5

# Numerical calculation of data table
# - varying p
cw_ModelB <- data.frame(p_vector)
colnames(cw_ModelB)[1] <- "p"
# - now calculating everything else over the vector of p
cw_ModelB$c <- c_B(cw_ModelB$p)
cw_ModelB$w <- w_B(cw_ModelB$p, b)
cw_ModelB$Ntotal <- N_eq14(cw_ModelB$p, b, total_area, AAw, AAc, aw, ac)
cw_ModelB$yield <- yield_B(cw_ModelB$p, b, total_area, AAw, AAc, aw, ac)
cw_ModelB$Nattrby_c <- N_eq14_c(cw_ModelB$p, b, total_area, AAw, AAc, aw, ac)
cw_ModelB$Nattrby_w <- N_eq14_w(cw_ModelB$p, b, total_area, AAw, AAc, aw, ac)

# Overall fig settings -------------
par(mfrow=c(1,2))
# bottom, left, top, right
par(mar=c(0.5, 3, 0.5, 0.5), oma=c(4, 0, 0, 0), mgp=c(3, 1, 0), las=1)
# Two graphs with differing attractiveness of w
# Fig. 3a Fixed area, identical resources ------------------
# Set up the frame
plot(NULL 
     , ylim=c(0,Flim_B)
     , xlim=c(0,1)
     , pch=17
     , cex=pointsize
     , xlab=""
     , ylab=""
)
# How bees divide up
# Fw
polygon(c(cw_ModelB$p, rev(cw_ModelB$p)), c(cw_ModelB$yield, rev(cw_ModelB$Ntotal))
#        , data = cw_ModelB
        , col = colorlist[5]
        , lty = 0
)
# Fc
polygon(c(cw_ModelB$p, rev(cw_ModelB$p)), c(c(1:resolution)*0, rev(cw_ModelB$yield))
#        , data = cw_ModelB
        , col = colorlist[2]
        , lty = 0
)
points(yield ~ p
       , data = cw_ModelB
       , col=colorlist[2]
       , pch=19
       , cex=pointsize)
# Where bees come from
points(Nattrby_w ~ p
       , data = cw_ModelB
       , col=colorlist[6]
       , pch=1
       , cex=pointsize
       )
points(Nattrby_c ~ p
       , data = cw_ModelB
       , col=colorlist[3]
       , pch=1
       , cex=pointsize
       )
# Total bees
points(Ntotal ~ p
       , data = cw_ModelB
       , col=colorlist[1]
       , pch=19
       , cex=pointsize)
legend("topright"
       , legend = c("Total area * bees (N)", "'Yield' for c (f)", "Bees attracted by c", "'Yield' for w", "Bees attracted by w")
       , col = c(colorlist[1], colorlist[2], colorlist[3], colorlist[5], colorlist[6])
       , pch = 19
       #, title= paste("Parameters b=", b, ", aw=", aw, sep="")
       , bty = "n")
# Overall x axis label
mtext(text = "Proportion of area that is w (p)",
      side = 1, line = 2, cex = 1, las = 1, outer = TRUE)
# Plotting horizontal and vertical lines at optimal p and yield
optimalY <- max(yield_B(p_vector, b, total_area, AAw, AAc, aw, ac))
optimalp <- p_vector[which(yield_B(p_vector, b, total_area, AAw, AAc, aw, ac)==optimalY)]
abline(h=optimalY, lty=2, lwd=2, col="grey")
abline(v=optimalp, lty=2, lwd=2, col="grey")
# Additional labels
text(0.1, optimalY*1.1, adj=0, paste("Optimal p*=", round(optimalp, 2), " \ngives f*=", round(optimalY,1), sep=""))
text(0.1, 0.98*Flim_B, "a", cex=2, adj=c(0,0))


# Parameters panel b ----------------------------
AAw <- 100
AAc <- 1
Flim_B <- 100

# Numerical calculation of data table
cw_ModelB <- data.frame(p_vector)
colnames(cw_ModelB)[1] <- "p"
cw_ModelB$c <- c_B(cw_ModelB$p)
cw_ModelB$w <- w_B(cw_ModelB$p, b)
cw_ModelB$Ntotal <- N_eq14(cw_ModelB$p, b, total_area, AAw, AAc, aw, ac)
cw_ModelB$yield <- yield_B(cw_ModelB$p, b, total_area, AAw, AAc, aw, ac)
cw_ModelB$Nattrby_c <- N_eq14_c(cw_ModelB$p, b, total_area, AAw, AAc, aw, ac)
cw_ModelB$Nattrby_w <- N_eq14_w(cw_ModelB$p, b, total_area, AAw, AAc, aw, ac)
# Fig. 3b Fixed area, one resource more attractive ----------------------------
plot(NULL
     , ylim=c(0,Flim_B)
     , xlim=c(0,1)
     , pch=17
     , cex=pointsize
)
# How bees divide up
# Fw
polygon(c(cw_ModelB$p, rev(cw_ModelB$p)), c(cw_ModelB$yield, rev(cw_ModelB$Ntotal))
        #        , data = cw_ModelB
        , col = colorlist[5]
        , lty = 0
)
# Fc
polygon(c(cw_ModelB$p, rev(cw_ModelB$p)), c(c(1:resolution)*0, rev(cw_ModelB$yield))
        #        , data = cw_ModelB
        , col = colorlist[2]
        , lty = 0
)
points(yield ~ p
       , data = cw_ModelB
       , col=colorlist[2]
       , pch=19
       , cex=pointsize)
# Where bees come from
points(Nattrby_w ~ p
       , data = cw_ModelB
       , col=colorlist[6] #alpha 0.01 to make it light color
       , pch=1
       , cex=pointsize
)
points(Nattrby_c ~ p
       , data = cw_ModelB
       , col=colorlist[3]
       , pch=1
       , cex=pointsize
)
# Total bees
points(Ntotal ~ p
       , data = cw_ModelB
       , col=colorlist[1]
       , pch=19
       , cex=pointsize
       )
# Additional labels and lines
optimalY <- max(yield_B(p_vector, b, total_area, AAw, AAc, aw, ac))
optimalp <- p_vector[which(yield_B(p_vector, b, total_area, AAw, AAc, aw, ac)==optimalY)]
abline(h=optimalY, lty=2, lwd=2, col="grey")
abline(v=optimalp, lty=2, lwd=2, col="grey")
text(optimalp*1.1, optimalY*1.25, adj=0, paste("Optimal p*=", round(optimalp, 2), " \ngives f*=", round(optimalY,1), sep=""))
text(0, 0.98*Flim_B, "b", cex=2, adj=c(0,0))

  




