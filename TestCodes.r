rm(list=ls())

# 载入包
library(igraph)
library(ggplot2)
library(dplyr)

# 定义核心函数：Cusp 动力学
stoch_cusp <- function(N, x, b, a, s_O, maxwell_convention=FALSE) {
  dt <- 0.01
  x <- x - dt * (x^3 - b * x - a) + rnorm(N, 0, s_O)
  if (maxwell_convention) x[x * a < 0] <- -x[x * a < 0]
  x[is.nan(x)] <- 0
  return(x)
}

# 参数设定
delta_list <- c(0.1, 0.4, 0.7, 0.9)
N <- 400
Ni <- 15000
min_attention <- -0.5
attention_star <- 1
s_O <- 0.01
persuasion <- 2
r_min <- 0.1
deffuant_c <- Inf
info_update <- TRUE
sd_noise_information <- 0.0005

# 构图（SBM 网络）
make_network <- function(N, clusters=10, p_within=0.2, p_between=0.001) {
  pm <- matrix(p_between, clusters, clusters)
  diag(pm) <- p_within
  g <- sample_sbm(N, pref.matrix = pm, block.sizes = rep(N / clusters, clusters))
  return(g)
}

# 结果存储
attention_all <- data.frame()

# 主循环：每个 delta_attention 值
for (delta_attention in delta_list) {
  set.seed(123 + as.integer(delta_attention*100))

  # 初始化网络与参数
  g <- make_network(N)
  attention <- runif(N, 0, 0)
  information <- rnorm(N, 0.1, 0)
  opinion <- rnorm(N, 0, 0.01)
  
  for (i in 1:500) {
    opinion <- stoch_cusp(N, opinion, attention + min_attention, information, s_O)
  }

  # 邻居矩阵
  max_deg <- max(degree(g))
  m_neigh <- matrix(NA, N, max_deg)
  for (i in 1:N) {
    nb <- neighbors(g, i)
    len <- length(nb)
    m_neigh[i, 1:len] <- nb
  }

  # 交互过程
  for (iter in 1:Ni) {
    agent <- sample(1:N, 1, prob = attention + 0.01)  # 避免全为0
    neigh <- na.omit(m_neigh[agent, ])
    if (length(neigh) == 0) next
    partner <- sample(neigh, 1)

    # bounded confidence 条件
    if (abs(opinion[agent] - opinion[partner]) < deffuant_c) {
      # 信息更新
      r1 <- r_min + (1 - r_min) / (1 + exp(-persuasion * (attention[agent] - attention[partner])))
      r2 <- r_min + (1 - r_min) / (1 + exp(-persuasion * (attention[partner] - attention[agent])))

      information[agent] <- if (info_update) r1 * information[agent] + (1 - r1) * information[partner] else information[agent]
      information[partner] <- if (info_update) r2 * information[partner] + (1 - r2) * information[agent] else information[partner]

      # attention 增加
      attention[agent] <- attention[agent] + delta_attention * (2 * attention_star - attention[agent])
      attention[partner] <- attention[partner] + delta_attention * (2 * attention_star - attention[partner])
    }

    # attention 衰减 + 信息扰动
    attention <- attention - 2 * delta_attention * attention / N
    information <- information + rnorm(N, 0, sd_noise_information)

    # opinion 更新
    opinion <- stoch_cusp(N, opinion, attention + min_attention, information, s_O)
  }

  # 存储结果
  attention_all <- rbind(attention_all,
                         data.frame(delta_A = as.factor(delta_attention), A = attention))
}

# Violin plot 展示结果
ggplot(attention_all, aes(x = delta_A, y = A, fill = delta_A)) +
  geom_violin(trim = FALSE, color = "black") +
  theme_minimal(base_size = 14) +
  labs(title = "Violin Plot of Attention A under different ΔA",
       x = expression(delta[A]), y = "Attention (A)") +
  scale_fill_brewer(palette = "Pastel2")
