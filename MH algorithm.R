## ---- 讀檔並建 bigram 模型：alphabet = a..z + space (共 27) ----
train_bigram_from_text <- function(file, alpha = 0.5) {
  # 讀取整本書
  txt <- tolower(paste(readLines(file, warn = FALSE, encoding = "UTF-8"),
                       collapse = " "))
  # 只保留 a-z 與空白；其他符號都轉空白，並合併連續空白
  txt <- gsub("[^a-z ]", " ", txt)
  txt <- gsub(" +", " ", txt)
  txt <- trimws(txt)
  
  # 字元 → 索引：a=1,...,z=26, space=27
  map_idx <- function(ch) if (ch == " ") 27L else (utf8ToInt(ch) - utf8ToInt("a") + 1L)
  chs <- strsplit(txt, "", fixed = TRUE)[[1]]
  x <- vapply(chs, map_idx, integer(1))
  
  # 統計 bigram 次數 N_ij
  K <- 27L
  N <- matrix(0L, K, K)
  if (length(x) >= 2) {
    for (t in 2:length(x)) N[x[t-1], x[t]] <- N[x[t-1], x[t]] + 1L
  }
  
  # Laplace (加法) 平滑：N + alpha
  P <- sweep(N + alpha, 1, rowSums(N + alpha), "/")
  M <- log(P)  # 對數轉移機率
  list(P = P, M = M, N = N)
}

bg <- train_bigram_from_text("C:/Users/user/Desktop/David_Copperfield.txt", alpha = 0.5)
M_lang <- bg$M 

## ---- 密文 ----
cipher_text <- "abvunr dqr obveybdexrhgesbqhgfe xoseznbwqgpe xoseosevb ervegrsyegtgnaosgeokeybdesbqhgfeo 
ewyejno ovueybdnebjveabfgeybdepds exrhgewggvekrpoqornejo xe xgepr gnorqe rdux eove xosejggleeexbzge xr 
eybdegvcbye xgesdppgnearpzervfe xr e xgepg xbfseybdeqgrnvgfexgngejoqqewgexgqzkdqe beybdefbevb exgso r ge beabv 
ra epgeokeybdexrhgervyemdgs obveove xgekd dng"

## ---- 基本映射與 bigram 計數（密文端） ----
alphabet <- c(letters, " ")
idx_of <- function(ch) if (ch == " ") 27L else match(ch, letters)
c_idx <- vapply(strsplit(cipher_text, "")[[1]], idx_of, integer(1))

K <- 27L
N_c <- matrix(0L, K, K)
if (length(c_idx) >= 2) {
  for (t in 2:length(c_idx)) N_c[c_idx[t-1], c_idx[t]] <- N_c[c_idx[t-1], c_idx[t]] + 1L
}

## ---- 對數似然（給定排列 perm）----
## perm 的意義：把「密文字母 i」映成「語言模型字母 perm[i]」
ll_of_perm <- function(perm, N_c, M_lang) {
  M_perm <- M_lang[perm, perm, drop = FALSE]
  sum(N_c * M_perm)
}

## ---- 將目前最佳排列解碼為可讀文字（預覽前 n 字） ----
decode_preview <- function(perm, c_idx, n = 120) {
  inv_map <- c(letters, " ")[perm]        # 密文 i → 明文字元
  s <- paste(inv_map[c_idx], collapse = "")
  substr(s, 1, min(nchar(s), n))
}

## ---- Metropolis–Hastings 在排列空間 ----
perm <- sample.int(K)               # 初始化：隨機排列
cur_ll <- ll_of_perm(perm, N_c, M_lang)
best_perm <- perm; best_ll <- cur_ll
acc <- 0L

Niter <- 80000
print_every <- 100

for (t in 1:Niter) {
  # 對稱提議：交換兩個符號
  a <- b <- 1L
  while (a == b) { a <- sample.int(K, 1); b <- sample.int(K, 1) }
  prop <- perm; prop[c(a,b)] <- prop[c(b,a)]
  
  d <- ll_of_perm(prop, N_c, M_lang) - cur_ll
  if (log(runif(1)) < d) {
    perm <- prop; cur_ll <- cur_ll + d; acc <- acc + 1L
    if (cur_ll > best_ll) { best_ll <- cur_ll; best_perm <- perm }
  }
  if (t %% print_every == 0L) {
    cat(sprintf("iter=%5d | ll=%.2f | best=%.2f | acc=%.3f | preview: %s\n",
                t, cur_ll, best_ll, acc/t,
                decode_preview(best_perm, c_idx, 80)))
  }
}

cat("\n=== 結果 ===\n")
cat(sprintf("接受率: %.3f\n最佳對數似然: %.2f\n", acc/Niter, best_ll))
cat("最佳解碼預覽：\n", decode_preview(best_perm, c_idx, 200), "\n", sep = "")

