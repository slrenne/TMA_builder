set.seed(1)
# spatial transcriptomics
# the chosen design is 0.6 cores 31x14 matrix
# there are 10 cases 
cases <- LETTERS[1:10]
# there are 3 core types [color on the slide]
# red: tumor center (R) [color Red]
# blue: non tumor (B) [color Blue]
# white: tumor margin (W) [color white]
core_type <- c('R','B','W')
# TMA_cores <- expand.grid(cases,core_type)
# TMA_cores <- paste(TMA_cores[,1],TMA_cores[,2],sep = '_')
# core_prob <- c(rep(0.05, times = 10), rep(0.01, times = 10),rep(0.04, times = 10))

# we will use a 31x14 TMA design
TMA_map <- matrix(nrow = 14, ncol = 31) #prepare a matrix to hold the TMA structure
# 'E' will be the empty columns and row to orient the TMA
TMA_map[,c(16,27)] <- 'E' 
TMA_map[11,] <- 'E'
TMA_map
# sum(is.na(TMA_map)) # 377
to_fill <- vector(length = sum(is.na(TMA_map)))

# ensure that there are at least 3 non neoplastic cores
to_fill[1:30] <- rep(paste(cases,core_type[2],sep='_'), times = 3) 

# ensure that there are at least 15 margin neoplastic cores
to_fill[31:210] <- rep(paste(cases,core_type[3],sep='_'), times = 18) 

# randomly fill the rest with tumor center
for(i in 211:length(to_fill)) { to_fill[i] <- 
  sample(paste(cases,core_type[1],sep='_'),
         size = 1, 
         replace = FALSE)}

pdf('samples.pdf', width = 9, height = 9)
plot(NULL, xlim = c(0,31), ylim = c(0,23), 
     xaxt = "n", yaxt = "n", xlab = "Samples", ylab = "nuber of cores")
lines(table(to_fill), lwd = 5, col = 1:3)
legend('topleft', lwd = 5, col = 1:3,  bty = "n",
       legend = c("Non tumor", "Tumor center", "Tumor Margin"))
axis(side = 1, at = 1:10*3 -1, labels = LETTERS[1:10])
axis(side = 2, at = 1:20, labels = 1:20)
dev.off()

# test the length of the vector and the empty spaces
sum(is.na(TMA_map)) == length(to_fill) # 377

to_fill_r <- sample(to_fill, 
                    size = length(to_fill),
                    replace = FALSE)

idx <- which(is.na(TMA_map))
TMA_v <- vector(length = 31*14)
for( i in 1:length(to_fill_r)) TMA_v[idx[i]] <- to_fill_r[i]
TMA <- matrix(TMA_v, nrow = 14, ncol = 31)
TMA[,c(16,27)] <- 'E' 
TMA[11,] <- 'E'


write.csv(TMA, 'map.csv')
pdf('map.pdf', width = 16, height = 11)
plot(expand.grid(1:31,1:14), cex = 6, 
     xlab='', ylab='', xaxt = 'n', yaxt = 'n')
axis(c(1),1:31)
axis(c(3),1:31)
axis(c(2),1:14,14:1)
axis(c(4),1:14,14:1)
for(c in 1:31){
  for(r in 1:14){
    text(x = c,
         y = 15-r,
         labels = TMA[r,c])}}
dev.off()

# Changes made during TMA building
TMA_t <- TMA
TMA_t[7,6] <- "H_R"
TMA_t[7,23] <- "H_R"
TMA_t[7,25] <- "E_R"
TMA_t[7,29] <- "E_R"
TMA_t[7,30] <- "E_R"

TMA_t[8,7] <- "E_R"
TMA_t[8,8] <- "H_R"
TMA_t[8,24] <- "H_R"

TMA_t[9,10] <- "H_R"
TMA_t[9,25] <- "H_R"
TMA_t[9,30] <- "I_R"

TMA_t[10,5] <- "H_R"
TMA_t[10,12] <- "I_R"
TMA_t[10,19] <- "I_R"
TMA_t[10,25] <- "H_R"

TMA_t[12,1] <- "I_R"
TMA_t[12,6] <- "I_R"
TMA_t[12,29] <- "J_R"
TMA_t[12,30] <- "J_B"

TMA_t[13,2] <- "H_R"
TMA_t[13,13] <- "I_R"
TMA_t[13,14] <- "I_B"
TMA_t[13,21] <- "H_R"
TMA_t[13,25] <- "J_B"
TMA_t[13,29] <- "E_R"
TMA_t[13,30] <- "H_W"

TMA_t[14,2] <- "H_B"
TMA_t[14,3] <- "E_R"
TMA_t[14,8] <- "J_B"
TMA_t[14,10] <- "E_W"
TMA_t[14,12] <- "E"
TMA_t[14,14] <- "J_B"
TMA_t[14,15] <- "E_R"
TMA_t[14,23] <- "J_R"
TMA_t[14,25] <- "E"

pdf('map_t.pdf', width = 16, height = 11)
plot(expand.grid(1:31,1:14), cex = 6, 
     xlab='', ylab='', xaxt = 'n', yaxt = 'n')
axis(c(1),1:31)
axis(c(3),1:31)
axis(c(2),1:14,14:1)
axis(c(4),1:14,14:1)
for(c in 1:31){
  for(r in 1:14){
    text(x = c,
         y = 15-r,
         labels = TMA_t[r,c])}}
dev.off()

# scanning with geomx
# column missing 30, 31
# row missing 14

TMA_g <- TMA_t[-13,-c(30,31)]
to_lines <- c(table(TMA_g)[1:9],
           0, # cutting out the missing column and the missing rows all 3 D_B are absent
           table(TMA_g)[10:11], #[12] represents empty spots
           table(TMA_g)[13:30])
attr(to_lines,"names")[10] <- "D_B"
pdf('samples_g.pdf', width = 9, height = 9)
plot(NULL, xlim = c(0,31), ylim = c(0,19), 
     xaxt = "n", yaxt = "n", xlab = "Samples", ylab = "nuber of cores")
segments(x0=1:30, y0=0, y1=to_lines, lwd = 5, col = 1:3)
legend('topright', lwd = 5, col = 1:3,  bty = "n",
       legend = c("Non tumor", "Tumor center", "Tumor Margin"))
axis(side = 1, at = 1:10*3 -1, labels = LETTERS[1:10])
axis(side = 2, at = 1:20, labels = 1:20)
dev.off()

TMA_g <- TMA_t[-13,-9:-31]
to_lines <- c(table(TMA_g)[1:6],
              0, # cutting out the missing column and the missing rows all 3 D_B are absent
              table(TMA_g)[7:8],
              0,
              table(TMA_g)[9:10],
              0,
              table(TMA_g)[12:14],
              0,#[12] represents empty spots
              table(TMA_g)[15:21],
              0,#[12] represents empty spots
              table(TMA_g)[22:26])
attr(to_lines,"names")[10] <- "D_B"
pdf('samples_g.pdf', width = 9, height = 9)

plot(NULL, xlim = c(0,31), ylim = c(0,19), 
     xaxt = "n", yaxt = "n", xlab = "Samples", ylab = "nuber of cores")
segments(x0=1:30, y0=0, y1=to_lines, lwd = 5, col = 1:3)
legend('topright', lwd = 5, col = 1:3,  bty = "n",
       legend = c("Non tumor", "Tumor center", "Tumor Margin"))
axis(side = 1, at = 1:10*3 -1, labels = LETTERS[1:10])
axis(side = 2, at = 1:20, labels = 1:20)
dev.off()
