require(ggplot2)
require(gridExtra)



# Table III

# heart
mass.heart.gm <- c(0.73, 51, 80, 120, 290, 320, 320, 330, 330, 335, 
                   370, 375, 2665, 3000, 2805, 2810, 2857, 2920)
mb.heart.pergm <- c(0.91, 2.1, 1.7, 1.7, NA, NA, 1.4, 1.1, 1.1, 1.0, 
                    NA, 1.3, 4.0, 4.7, 5.1, 4.7, 4.5, 4.4)
cytc.heart.pergm <- c(0.447, 0.230, 0.281, 0.228, 0.108, 0.148, 0.128,
                      0.137, 0.127, 0.146, 0.158, 0.138, 0.148, 0.212,
                      0.152, 0.182, 0.171, 0.187)

# skeletal muscle
mass.skeletal.gm <- c(0.113, 2.286, 3.557, 5.328, 26.04, 28.14, 28.56,
                      29.40, 29.40, 30.24, 31.92, 32.76, 200, 205, 207,
                      216, 225, 225)*(10^3)
mb.skeletal.pergm <- c(0.89, 5.1, 2.1, 2.8, NA, NA, 1.4, 1.1, 1.0, 1.0,
                       NA, 1.4, 6.2, 6.6, 8.2, 7.8, 7.0, 8.3)
cytc.skeletal.pergm <- c(0.098, 0.048, 0.056, 0.052, 0.018, 0.024, 0.020,
                         0.023, 0.018, 0.022, 0.028, 0.024, 0.046, 0.095,
                         0.042, 0.057, 0.071, 0.059)

# species
species <- c("R", "D", "D", "D", "M", "M", 
             "M", "M", "M", "M", "M", "M", 
             "H", "H", "H", "H", "H", "H")

d <- data.frame(mass.heart.gm = mass.heart.gm,
                mb.heart.pergm = mb.heart.pergm,
                cytc.heart.pergm = cytc.heart.pergm,
                mb.heart.mg = mass.heart.gm*mb.heart.pergm,
                cytc.heart.mg = mass.heart.gm*cytc.heart.pergm,
                mass.skeletal.gm = mass.skeletal.gm,
                mb.skeletal.pergm = mb.skeletal.pergm,
                cytc.skeletal.pergm = cytc.skeletal.pergm,
                mb.skeletal.mg = mass.skeletal.gm*mb.skeletal.pergm,
                cytc.skeletal.mg = mass.skeletal.gm*cytc.skeletal.pergm,
                species = species)



grid.arrange(
  ggplot(data = d)+
    geom_point(aes(x = mass.heart.gm, y = mb.heart.pergm, col = species))+
    theme_bw(base_size = 9),
  ggplot(data = d)+
    geom_point(aes(x = mass.heart.gm, y = mb.heart.mg, col = species))+
    theme_bw(base_size = 9)+
    scale_x_log10()+
    scale_y_log10()+
    geom_abline(slope = 1.21, intercept = -0.228),
  ggplot(data = d)+
    geom_point(aes(x = mass.heart.gm, y = cytc.heart.mg, col = species))+
    theme_bw(base_size = 9)+
    scale_x_log10()+
    scale_y_log10()+
    geom_abline(slope = 0.899, intercept = -0.5),
  ggplot(data = d)+
    geom_point(aes(x = mass.skeletal.gm, y = mb.skeletal.pergm, col = species))+
    theme_bw(base_size = 9),
  ggplot(data = d)+
    geom_point(aes(x = mass.skeletal.gm, y = mb.skeletal.mg, col = species))+
    theme_bw(base_size = 9)+
    scale_x_log10()+
    scale_y_log10()+
    geom_abline(slope = 1.245, intercept = -0.633),
  ggplot(data = d)+
    geom_point(aes(x = mass.skeletal.gm, y = cytc.skeletal.mg, col = species))+
    theme_bw(base_size = 9)+
    scale_x_log10()+
    scale_y_log10()+
    geom_abline(slope = 0.972, intercept = -1.294),
  nrow = 2)


# heart
summary(lm(log10(d$mb.heart.mg)~log10(d$mass.heart.gm)))
summary(lm(log10(d$cytc.heart.mg)~log10(d$mass.heart.gm)))

# skeletal muscle
summary(lm(log10(d$mb.skeletal.mg)~log10(d$mass.skeletal.gm)))
summary(lm(log10(d$cytc.skeletal.mg)~log10(d$mass.skeletal.gm)))







# Table IV
species <- c("R", "D", "D", "M", "He", "Ho", "Ho")
mass.gm <- c(0.250, 9.88, 6.35, 70, 182, 500, 455)*(10^3)
skeletal.mass.gm <- c(0.113, 3.56, 2.29, 29.40, 81.9, 225, 205)*(10^3)

hb.gm <- c(3.19, 138.3, 88.9, 912.8, 2215, 5805, NA)
mb.gm <- c(0.101, 7.5, 11.7, 34.7, 307, 1867.5, 1345.2)
cytc.gm <- c(0.0144, 0.249, 0.137, 0.781, 1.24, 16.6, 24.3)


t <- data.frame(species = species,
                mass.gm = mass.gm,
                skeletal.mass.gm = skeletal.mass.gm,
                hb.gm = hb.gm,
                mb.gm = mb.gm,
                cytc.gm = cytc.gm)

summary(lm(log10(t$hb.gm)~log10(t$mass.gm)))
summary(lm(log10(t$cytc.gm)~log10(t$mass.gm)))
summary(lm(log10(t$mb.gm)~log10(t$mass.gm)))
rm(species, mass.gm, skeletal.mass.gm, hb.gm, mb.gm, cytc.gm)

grid.arrange(
  ggplot(data = t)+
    geom_smooth(aes(x = mass.gm, y = hb.gm), method = 'lm')+
    geom_point(aes(x = mass.gm, y = hb.gm, col = species))+
    theme_bw(base_size = 9)+
    theme(legend.position = "top")+
    scale_x_log10()+
    scale_y_log10()+
    geom_abline(slope = 1, intercept = -1.825, col = "darkgray", linetype = "dashed"),
  ggplot(data = t)+
    geom_smooth(aes(x = mass.gm, y = mb.gm), method = 'lm')+
    geom_point(aes(x = mass.gm, y = mb.gm, col = species))+
    theme_bw(base_size = 9)+
    theme(legend.position = "top")+
    scale_x_log10()+
    scale_y_log10()+
    geom_abline(slope = 1, intercept = -3.937, col = "darkgray", linetype = "dashed"),
  ggplot(data = t)+
    geom_smooth(aes(x = mass.gm, y = cytc.gm), method = 'lm')+
    geom_point(aes(x = mass.gm, y = cytc.gm, col = species))+
    theme_bw(base_size = 9)+
    theme(legend.position = "top")+
    scale_x_log10()+
    scale_y_log10()+
    geom_abline(slope = 1, intercept = -4.239, col = "darkgray", linetype = "dashed"),
  nrow = 1)
