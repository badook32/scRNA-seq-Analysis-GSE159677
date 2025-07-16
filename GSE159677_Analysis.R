# 1 준비 (Preparation)

# 1.1 R 패키지 로드
library(tidyverse)
library(Seurat)
library(scran)
library(patchwork)
library(viridis)
library(ggforce)
library(gghalves)
library(ggridges)
library(scDblFinder)
library(SingleR)
library(cerebroApp)


# 1.2 색상 정의
custom_colors <- list()

# 네덜란드 스타일 색상 팔레트
colors_dutch <- c(
  '#FFC312','#C4E538','#12CBC4','#FDA7DF','#ED4C67',
  '#F79F1F','#A3CB38','#1289A7','#D980FA','#B53471',
  '#EE5A24','#009432','#0652DD','#9980FA','#833471',
  '#EA2027','#006266','#1B1464','#5758BB','#6F1E51'
)

# 스페인 스타일 색상 팔레트
colors_spanish <- c(
  '#40407a','#706fd3','#f7f1e3','#34ace0','#33d9b2',
  '#2c2c54','#474787','#aaa69d','#227093','#218c74',
  '#ff5252','#ff793f','#d1ccc0','#ffb142','#ffda79',
  '#b33939','#cd6133','#84817a','#cc8e35','#ccae62'
)

# 이산형 데이터(범주형 데이터)용 색상 지정
custom_colors$discrete <- c(colors_dutch, colors_spanish)

# 세포 주기 단계별 색상 설정
custom_colors$cell_cycle <- setNames(
  c('#45aaf2', '#f1c40f', '#e74c3c', '#7f8c8d'),
  c('G1',      'S',       'G2M',     '-')
)

# 1.3 보조 함수 (Helper functions)

# 1.3.1 숫자 포맷팅 함수
niceFormat <- function(number) {
  formatC(number, format = 'f', big.mark = ',', digits = 0)
}

# 1.3.2 10x Genomics에서 생성된 전사체 카운트 행렬 로드 함수
load10xData <- function(
    sample_name, # 샘플명
    path, # 데이터 파일이 저장된 경로
    max_number_of_cells = NULL # 선택적: 최대 셀 개수 제한
) {
  
  ## 전사체 카운트 행렬 읽기
  transcript_counts <- list.files(
    path,
    recursive = TRUE,
    full.names = TRUE
  ) %>%
    .[grep(., pattern = 'mtx')] %>%
    Matrix::readMM()
  
  ## 세포 바코드 파일을 읽고, 이를 전사체 카운트 행렬의 열 이름으로 설정
  colnames(transcript_counts) <- list.files(
    path,
    recursive = TRUE,
    full.names = TRUE
  ) %>%
    .[grep(., pattern = 'barcodes')] %>%
    .[grep(., pattern = 'tsv')] %>%
    read_tsv(col_names = FALSE, col_types = cols()) %>%
    pull(1) %>%
    gsub(., pattern = '-[0-9]{1,2}', replacement = '') %>%
    paste0(., '-', sample_name)
  
  ## 유전자 이름 파일을 읽기
  gene_names <- list.files(
    path,
    recursive = TRUE,
    full.names = TRUE
  ) %>%
    .[grep(., pattern = 'genes|features')] %>%
    read_tsv(col_names = FALSE, col_types = cols()) %>%
    pull(ifelse(ncol(.) > 1, 2, 1))
  
  ## 동일 유전자의 여러 전사체 데이터를 합쳐야 하는 경우 처리
  if ( any(duplicated(gene_names)) == TRUE ) {
    
    ## 중복된 유전자 이름 목록 가져오기
    duplicated_gene_names <- unique(gene_names[which(duplicated(gene_names))])
    
    ## 로그 메시지 출력
    message(
      glue::glue(
        "{format(Sys.time(), '[%Y-%m-%d %H:%M:%S]')} 중복된 유전자 {length(duplicated_gene_names)}개의 카운트를 합산 중."
      )
    )
    
    ## 중복된 각 유전자에 대해 처리 수행
    for ( i in duplicated_gene_names ) {
      
      ## 현재 유전자에 대한 전사체 카운트 추출
      tmp_counts <- transcript_counts[which(gene_names == i),]
      
      ## 적어도 2개 이상의 행이 존재하는 경우만 처리
      if ( nrow(tmp_counts) >= 2 ) {
        
        ## 각 세포에서의 카운트 합산
        tmp_counts_sum <- Matrix::colSums(tmp_counts)
        
        ## 기존 전사체 행렬에서 원래 카운트 제거 및 유전자 이름 제거
        transcript_counts <- transcript_counts[-c(which(gene_names == i)),]
        gene_names <- gene_names[-c(which(gene_names == i))]
        
        ## 합산된 카운트를 전사체 행렬의 마지막 행에 추가
        ## 유전자 이름도 마지막에 추가
        transcript_counts <- rbind(transcript_counts, tmp_counts_sum)
        gene_names <- c(gene_names, i)
      }
    }
  }
  
  ## 유전자 이름을 전사체 행렬의 행 이름으로 지정
  rownames(transcript_counts) <- gene_names
  
  ## 만약 `max_number_of_cells`보다 많은 세포가 존재하면, 일부만 샘플링하여 유지
  if (
    !is.null(max_number_of_cells) &&
    is.numeric(max_number_of_cells) &&
    max_number_of_cells < ncol(transcript_counts)
  ) {
    temp_cells_to_keep <- base::sample(1:ncol(transcript_counts), max_number_of_cells)
    transcript_counts <- transcript_counts[,temp_cells_to_keep]
  }
  
  ##
  message(
    glue::glue(
      "{format(Sys.time(), '[%Y-%m-%d %H:%M:%S]')} 샘플 {sample_name}의 카운트 로드 완료: ",
      "{niceFormat(ncol(transcript_counts))}개의 세포와 ",
      "{niceFormat(nrow(transcript_counts))}개의 유전자."
    )
  )
  
  ##
  return(transcript_counts)
}

# 1.3.3 유전자 이름이 동일한지 확인하는 함수
sameGeneNames <- function(counts) {
  
  ## 각 카운트 행렬에서 유전자 이름을 저장할 리스트 생성
  gene_names <- list()
  
  ## 각 카운트 행렬을 순회
  for ( i in names(counts) ) {
    
    ## 유전자 이름 추출
    gene_names[[i]] <- rownames(counts[[i]])
  }
  
  ## 첫 번째 카운트 행렬을 기준으로 모든 유전자 이름이 동일한지 확인
  return(all(sapply(gene_names, FUN = identical, gene_names[[1]])))
}

# 1.3.4 샘플별 nCount 및 nFeature 통계 테이블 생성 함수
averagenCountnFeature <- function(cells) {
  output <- tibble(
    'sample' = character(), # 샘플명
    'mean(nCount)' = numeric(), # nCount 평균
    'median(nCount)' = numeric(), # nCount 중앙값
    'mean(nFeature)' = numeric(), # nFeature 평균
    'median(nFeature)' = numeric() # nFeature 중앙값
  )
  
  ## 각 샘플을 순회하면서 통계값 계산
  for ( i in levels(cells$sample) ) {
    tmp <- tibble(
      'sample' = i,
      'mean(nCount)' = cells %>% filter(sample == i) %>% pull(nCount) %>% mean(),
      'median(nCount)' = cells %>% filter(sample == i) %>% pull(nCount) %>% median(),
      'mean(nFeature)' = cells %>% filter(sample == i) %>% pull(nFeature) %>% mean(),
      'median(nFeature)' = cells %>% filter(sample == i) %>% pull(nFeature) %>% median()
    )
    output <- bind_rows(output, tmp)
  }
  
  return(output)
}

# 1.4 출력 디렉토리 생성
dir.create('data')  # 데이터 저장 폴더 생성
dir.create('plots') # 시각화 결과 저장 폴더 생성

# ----------------------------------------------

setwd("~/Desktop/bioinformatics/AS/GSE159677_RAW(scRNA,ACvsPA)/GSE159677")

# ----------------------------------------------

# 2 데이터 로드 (Load data)
# 2.1 샘플별 데이터 로드 (Sample-wise)
sample_names <- list.dirs('raw_data', recursive = FALSE) %>% basename()

# AC 샘플과 PA 샘플을 구분
AC_samples <- sample_names[grepl("AC", sample_names)]
PA_samples <- sample_names[grepl("PA", sample_names)]

# 전사체 카운트 데이터를 저장할 리스트 생성
transcripts <- list(
  raw = list()
)

# AC 샘플 데이터 로드
for (i in AC_samples) {
  transcripts$raw[[i]] <- load10xData(
    i, paste0('raw_data/', i, '/'),
    max_number_of_cells = 2000 # 최대 2000개의 세포만 유지
  )
}

# PA 샘플 데이터 로드
for (i in PA_samples) {
  transcripts$raw[[i]] <- load10xData(
    i, paste0('raw_data/', i, '/'),
    max_number_of_cells = 2000 # 최대 2000개의 세포만 유지
  )
}

# 2.2 샘플 병합 (Merge samples)
sameGeneNames(transcripts$raw) # 모든 샘플에서 유전자 이름이 동일한지 확인
# TRUE 반환 시 병합 가능

# 모든 샘플의 전사체 카운트 행렬을 병합 (AC 샘플 + PA 샘플)
transcripts$raw$merged <- do.call(cbind, transcripts$raw)

# 병합된 행렬의 차원 확인
dim(transcripts$raw$merged)

# 유전자명을 명시적으로 컬럼으로 변환
valid_samples <- transcripts$raw[names(transcripts$raw) != "merged"]  # 병합 데이터 제외
valid_samples <- valid_samples[!sapply(valid_samples, is.null)]  # NULL 값 제외

# 각 샘플에 대해 처리
for (i in seq_along(valid_samples)) {
  # "gene" 컬럼이 이미 존재하는 경우 삭제
  if ("gene" %in% colnames(valid_samples[[i]])) {
    valid_samples[[i]] <- valid_samples[[i]] %>% select(-gene)
  }
  
  # 유전자명을 "gene" 컬럼으로 변환
  valid_samples[[i]] <- valid_samples[[i]] %>%
    as.data.frame() %>%
    rownames_to_column(var = "gene")
}

# 유전자별 발현값 합산
for (i in seq_along(valid_samples)) {
  valid_samples[[i]] <- valid_samples[[i]] %>%
    group_by(gene) %>%
    summarise(across(where(is.numeric), ~ sum(.x, na.rm = TRUE))) %>%
    ungroup()
}

# 유전자 기준으로 병합
transcripts$raw$merged <- Reduce(
  function(x, y) full_join(x, y, by = "gene"),  # gene을 기준으로 병합
  valid_samples
)

# NA 값을 0으로 변환
transcripts$raw$merged[is.na(transcripts$raw$merged)] <- 0

# 유전자명을 rownames로 변환
transcripts$raw$merged <- transcripts$raw$merged %>%
  column_to_rownames(var = "gene")

# 유전자 유실 여부 확인
original_genes <- unique(unlist(lapply(valid_samples, function(x) x$gene)))  # 원본 모든 샘플의 유전자 목록
merged_genes <- rownames(transcripts$raw$merged)  # 병합 후 유전자 목록

missing_genes <- setdiff(original_genes, merged_genes)

if (length(missing_genes) == 0) {
  cat("모든 유전자가 정상적으로 병합되었습니다.\n")
} else {
  cat("⚠ 병합 과정에서", length(missing_genes), "개의 유전자가 유실됨!\n")
  print(missing_genes[1:20])  # 일부 유실된 유전자 출력
}

# 최종 데이터 구조 확인
str(transcripts$raw$merged)


# ----------------------------------

# 3 품질 관리 (Quality control)

# 3.1 세포 필터링 (Filter cells)
# 3.1.1 데이터 준비 (샘플명을 AC / PA로 변환)
cells <- tibble(
  cell = colnames(transcripts$raw$merged),  # 세포 ID
  sample = colnames(transcripts$raw$merged) %>%
    sub(".*-Patient [1-3] ", "", .) %>%  # 'Patient X ' 앞부분(세포 ID 포함) 제거
    sub(" scRNA-seq$", "", .),  # ' scRNA-seq' 부분 제거하여 'AC' 또는 'PA'만 남김
  nCount = colSums(transcripts$raw$merged),  # 총 전사체 수 계산
  nFeature = colSums(transcripts$raw$merged != 0)  # 발현된 유전자 수 계산
)

# sample을 factor로 변환하여 AC, PA 순서 유지
cells$sample <- factor(cells$sample, levels = c("AC", "PA"))

# 결과 확인
print(table(cells$sample))  # AC와 PA 개수 확인


# 데이터 프레임 출력 (샘플 6000개, 4개 변수)
DataFrame(cells)
# DataFrame with 6000 rows and 4 columns
#                    cell   sample    nCount  nFeature
#             <character> <factor> <numeric> <integer>
# 1    AGTAGTCAGCTAAACA-A        A     11915       557
# 2    GCAAACTGTGTGCCTG-A        A      4745      1018
# 3    TTGGCAAAGGCAAAGA-A        A      4205      1027
# 4    CTACACCTCAACGAAA-A        A     11593       633
# 5    AGCCTAACAAGGTTCT-A        A      4756      1298
# ...                 ...      ...       ...       ...
# 5996 CCACCTAGTTAAAGAC-E        E      4941      1442
# 5997 AGGGAGTTCCTGTAGA-E        E      2405       699
# 5998 CATTCGCAGGTACTCT-E        E      2262       888
# 5999 TGATTTCAGGCTCAGA-E        E      2150       733
# 6000 CAGCTGGGTTATCGGT-E        E     21167      4030

# 미토콘드리아 유전자 목록 불러오기
mitochondrial_genes_here <- read_tsv(
  system.file('extdata/genes_mt_hg_name.tsv.gz', package = 'cerebroApp'),
  col_names = FALSE
) %>%
  filter(X1 %in% rownames(transcripts$raw$merged)) %>% # 전사체 행렬에 존재하는 유전자만 필터링
  pull(X1) # 유전자 이름 추출

# 각 세포의 미토콘드리아 유전자 비율(percent_MT) 계산
# transcripts$raw$merged에서 'gene' 열을 제외한 행렬을 가져옵니다.
merged_data <- transcripts$raw$merged[, -1]  # gene 열 제외

# 1) merged_data와 cells의 세포 ID가 동일한지 확인
common_cells <- intersect(colnames(merged_data), cells$cell)

# 2) 공통 세포만 유지하여 정렬
merged_data <- merged_data[, common_cells]
cells <- cells %>% filter(cell %in% common_cells)

# 3) 미토콘드리아 유전자 목록 확인 후 존재하는 유전자만 선택
mitochondrial_genes <- rownames(merged_data)[grepl("^MT-", rownames(merged_data))]

# 4) 미토콘드리아 유전자의 발현 비율 계산
cells$percent_MT <- Matrix::colSums(merged_data[mitochondrial_genes, , drop = FALSE]) /
  Matrix::colSums(merged_data)

# 5) 결과 확인
print(dim(cells))  # [1] 11999  5 (12000 -> 11999로 일치)
print(head(cells$percent_MT))  # 정상적으로 값이 들어갔는지 확인

# 데이터 프레임 업데이트 (percent_MT 추가됨)
DataFrame(cells)
# DataFrame with 6000 rows and 5 columns
#                    cell   sample    nCount  nFeature percent_MT
#             <character> <factor> <numeric> <integer>  <numeric>
# 1    AGTAGTCAGCTAAACA-A        A     11915       557 0.00260176
# 2    GCAAACTGTGTGCCTG-A        A      4745      1018 0.02971549
# 3    TTGGCAAAGGCAAAGA-A        A      4205      1027 0.02354340
# 4    CTACACCTCAACGAAA-A        A     11593       633 0.00517554
# 5    AGCCTAACAAGGTTCT-A        A      4756      1298 0.06391926
# ...                 ...      ...       ...       ...        ...
# 5996 CCACCTAGTTAAAGAC-E        E      4941      1442  0.0333940
# 5997 AGGGAGTTCCTGTAGA-E        E      2405       699  0.0191268
# 5998 CATTCGCAGGTACTCT-E        E      2262       888  0.0539346
# 5999 TGATTTCAGGCTCAGA-E        E      2150       733  0.0232558
# 6000 CAGCTGGGTTATCGGT-E        E     21167      4030  0.0313223

# 3.1.2 Doublets

library(Matrix)

# 'merged' 데이터를 Matrix 형식으로 변환
merged_matrix <- as(transcripts$raw$merged, "Matrix")

# SingleCellExperiment 객체 생성
sce <- SingleCellExperiment(assays = list(counts = merged_matrix))

# Doublets 탐지
sce <- scDblFinder(sce)

# 1. cells와 sce의 행 수가 동일한지 확인
cat("cells 행 수:", nrow(cells), "\n")
cat("sce 행 수:", ncol(sce), "\n")

# 2. 만약 행 수가 다르면, 하나를 제거하여 일치시킴
if (nrow(cells) != ncol(sce)) {
  message("행 수 불일치! 마지막 세포를 제거하여 일치시킵니다.")
  sce <- sce[, 1:nrow(cells)]  # 마지막 하나를 제거하여 동일한 행 수로 맞추기
}

# 3. multiplet_class 값을 cells에 추가
cells$multiplet_class <- colData(sce)$scDblFinder.class

# 결과 확인
head(cells)

p <- cells %>%
  filter(multiplet_class != 'singlet') %>%
  group_by(sample) %>%
  summarize(count = n()) %>%
  ggplot(aes(x = sample, y = count, fill = sample)) +
  geom_col(color = 'black') +
  theme_bw() +
  scale_fill_manual(values = custom_colors$discrete) +
  scale_x_discrete(limits = rev(levels(cells$sample))) +
  scale_y_continuous(name = 'Number of doublets', labels = scales::comma) +
  theme(
    axis.title.y = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.minor.y = element_blank(),
    legend.position = 'none'
  ) +
  coord_flip()

ggsave(
  'plots/qc_number_of_doublets_by_sample.png',
  p, height = 3, width = 6
)

# 'cells' 데이터에서 PA와 AC 그룹만 추출
cells$sample_group <- ifelse(grepl("PA", cells$sample), "PA", "AC")

# 각 sample_group별 nCount(총 전사체 수) 시각화
p1 <- ggplot(cells, aes(x = sample_group, y = nCount, fill = sample_group)) +
  geom_violin(draw_quantiles = c(0.5), scale = 'area', trim = FALSE) + # 중위값 표시
  theme_bw() +
  scale_fill_manual(values = custom_colors$discrete) +
  scale_x_discrete(limits = rev(levels(cells$sample_group))) +
  scale_y_continuous(labels = scales::comma) +
  labs(title = 'Number of transcripts', subtitle = 'linear scale') +
  theme(
    axis.title = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.minor.y = element_blank(),
    legend.position = 'none'
  ) +
  coord_flip()

# 로그 스케일로 변환한 nCount 시각화
p2 <- ggplot(cells, aes(x = sample_group, y = nCount, fill = sample_group)) +
  geom_violin(draw_quantiles = c(0.5), scale = 'area', trim = FALSE) +
  theme_bw() +
  scale_fill_manual(values = custom_colors$discrete) +
  scale_x_discrete(limits = rev(levels(cells$sample_group))) +
  scale_y_log10(labels = scales::comma) + # 로그 스케일 적용
  labs(title = 'Number of transcripts', subtitle = 'log-scale') +
  theme(
    axis.title = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.minor.y = element_blank(),
    legend.position = 'none'
  ) +
  coord_flip()

# 각 sample_group별 nFeature(발현된 유전자 수) 시각화
p3 <- ggplot(cells, aes(x = sample_group, y = nFeature, fill = sample_group)) +
  geom_violin(draw_quantiles = c(0.5), scale = 'area', trim = FALSE) +
  theme_bw() +
  scale_fill_manual(values = custom_colors$discrete) +
  scale_x_discrete(limits = rev(levels(cells$sample_group))) +
  scale_y_continuous(labels = scales::comma) +
  labs(title = 'Number of expressed genes', subtitle = 'linear scale') +
  theme(
    axis.title = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.minor.y = element_blank(),
    legend.position = 'none'
  ) +
  coord_flip()

# 로그 스케일로 변환한 nFeature 시각화
p4 <- ggplot(cells, aes(x = sample_group, y = nFeature, fill = sample_group)) +
  geom_violin(draw_quantiles = c(0.5), scale = 'area', trim = FALSE) +
  theme_bw() +
  scale_fill_manual(values = custom_colors$discrete) +
  scale_x_discrete(limits = rev(levels(cells$sample_group))) +
  scale_y_log10(labels = scales::comma) + # 로그 스케일 적용
  labs(title = 'Number of expressed genes', subtitle = 'log-scale') +
  theme(
    axis.title = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.minor.y = element_blank(),
    legend.position = 'none'
  ) +
  coord_flip()

# 그래프 저장 (nCount 및 nFeature 분포)
ggsave(
  'plots/qc_ncount_nfeature_by_sample_group.png',
  p1 + p3 + p2 + p4 + plot_layout(ncol = 2),
  height = 7, width = 10
)

# 중복 세포(doublet) 제거 (싱글 세포만 유지)
cells <- filter(cells, multiplet_class == 'singlet')

# 3.1.3 필터링 전 데이터 시각화 (Before filtering)
# 옵션 1: 각 샘플별 전사체 수 및 발현된 유전자 수 분포 확인
# 각 샘플별 nCount, nFeature, percent_MT 시각화
# sample 컬럼을 PA와 AC로 나누기
cells$sample_group <- ifelse(grepl("PA", cells$sample), "PA", "AC")

# 1. 선형 스케일 전사체 수 (nCount)
p1 <- ggplot(cells, aes(x = sample_group, y = nCount, fill = sample_group)) +
  geom_violin(draw_quantiles = c(0.5), scale = 'area', trim = FALSE) + # 중앙값 표시
  theme_bw() +
  scale_fill_manual(values = custom_colors$discrete) +
  scale_x_discrete(limits = rev(levels(cells$sample_group))) +
  scale_y_continuous(labels = scales::comma) +
  labs(title = 'Number of transcripts', subtitle = 'linear scale') +
  theme(
    axis.title = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.minor.y = element_blank(),
    legend.position = 'none'
  ) +
  coord_flip()

# 2. 로그 스케일 전사체 수 (nCount)
p2 <- ggplot(cells, aes(x = sample_group, y = nCount, fill = sample_group)) +
  geom_violin(draw_quantiles = c(0.5), scale = 'area', trim = FALSE) +
  theme_bw() +
  scale_fill_manual(values = custom_colors$discrete) +
  scale_x_discrete(limits = rev(levels(cells$sample_group))) +
  scale_y_log10(labels = scales::comma) +
  labs(title = 'Number of transcripts', subtitle = 'log-scale') +
  theme(
    axis.title = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.minor.y = element_blank(),
    legend.position = 'none'
  ) +
  coord_flip()

# 3. 선형 스케일 발현된 유전자 수 (nFeature)
p3 <- ggplot(cells, aes(x = sample_group, y = nFeature, fill = sample_group)) +
  geom_violin(draw_quantiles = c(0.5), scale = 'area', trim = FALSE) +
  theme_bw() +
  scale_fill_manual(values = custom_colors$discrete) +
  scale_x_discrete(limits = rev(levels(cells$sample_group))) +
  scale_y_continuous(labels = scales::comma) +
  labs(title = 'Number of expressed genes', subtitle = 'linear scale') +
  theme(
    axis.title = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.minor.y = element_blank(),
    legend.position = 'none'
  ) +
  coord_flip()

# 4. 로그 스케일 발현된 유전자 수 (nFeature)
p4 <- ggplot(cells, aes(x = sample_group, y = nFeature, fill = sample_group)) +
  geom_violin(draw_quantiles = c(0.5), scale = 'area', trim = FALSE) +
  theme_bw() +
  scale_fill_manual(values = custom_colors$discrete) +
  scale_x_discrete(limits = rev(levels(cells$sample_group))) +
  scale_y_log10(labels = scales::comma) +
  labs(title = 'Number of expressed genes', subtitle = 'log-scale') +
  theme(
    axis.title = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.minor.y = element_blank(),
    legend.position = 'none'
  ) +
  coord_flip()

# 5. 선형 스케일 미토콘드리아 전사체 비율 (percent_MT)
p5 <- ggplot(cells, aes(x = sample_group, y = percent_MT, fill = sample_group)) +
  geom_violin(draw_quantiles = c(0.5), scale = 'area', trim = FALSE) +
  theme_bw() +
  scale_fill_manual(values = custom_colors$discrete) +
  scale_x_discrete(limits = rev(levels(cells$sample_group))) +
  scale_y_continuous(labels = scales::percent) +
  labs(title = 'Percent MT transcripts', subtitle = 'linear scale') +
  theme(
    axis.title = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.minor.y = element_blank(),
    legend.position = 'none'
  ) +
  coord_flip()

# 필터링 전 nCount, nFeature, percent_MT의 분포 그래프 저장
ggsave(
  'plots/qc_histogram_ncount_nfeature_percentmt_before_filtering_1.png',
  p1 + p3 + p5 +
    p2 + p4 + plot_layout(ncol = 3),
  height = 7, width = 10
)

# 3.1.4 임곗값 정의 (Define thresholds)
# 각 샘플을 PA와 AC로 나누기
cells$sample_group <- ifelse(grepl("PA", cells$sample), "PA", "AC")

# 1. nCount(전사체 수) 임곗값 계산
median_nCount <- median(cells$nCount) # 중앙값: 4255
mad_nCount <- mad(cells$nCount) # MAD: 2630.874

# 2. nFeature(발현된 유전자 수) 임곗값 계산
median_nFeature <- median(cells$nFeature) # 중앙값: 803
mad_nFeature <- mad(cells$nFeature) # MAD: 544.1142

# 3. 미토콘드리아 전사체 비율(percent_MT) 임곗값 계산
median_percent_MT <- median(cells$percent_MT) # 중앙값: 0.02760973
mad_percent_MT <- mad(cells$percent_MT) # MAD: 0.02103674

# 임곗값 설정 (중앙값 + 5*MAD)
thresholds_nCount <- c(0, median_nCount + 5 * mad_nCount) # [0, 17409.37]
thresholds_nFeature <- c(0, median_nFeature + 5 * mad_nFeature) # [0, 3523.571]
thresholds_percent_MT <- c(0, median_percent_MT + 5 * mad_percent_MT) # [0, 0.1327934]

# 필터링 전 전사체 수 (nCount) 시각화 (선형 스케일)
p1 <- ggplot(cells, aes(x = sample_group, y = nCount, fill = sample_group)) +
  geom_violin(draw_quantiles = c(0.5), scale = 'area', trim = FALSE) +
  geom_hline(yintercept = median_nCount, color = 'black') +
  geom_hline(yintercept = thresholds_nCount[2], color = 'red') +
  theme_bw() +
  scale_fill_manual(values = custom_colors$discrete) +
  scale_x_discrete(limits = rev(levels(cells$sample_group))) +
  scale_y_continuous(labels = scales::comma) +
  labs(title = 'Number of transcripts', subtitle = 'linear scale') +
  theme(axis.title = element_blank(), legend.position = 'none') +
  coord_flip()

# 필터링 전 전사체 수 (nCount) 시각화 (로그 스케일)
p2 <- ggplot(cells, aes(x = sample_group, y = nCount, fill = sample_group)) +
  geom_violin(draw_quantiles = c(0.5), scale = 'area', trim = FALSE) +
  geom_hline(yintercept = median_nCount, color = 'black') +
  geom_hline(yintercept = thresholds_nCount[2], color = 'red') +
  theme_bw() +
  scale_fill_manual(values = custom_colors$discrete) +
  scale_x_discrete(limits = rev(levels(cells$sample_group))) +
  scale_y_log10(labels = scales::comma) +
  labs(title = 'Number of transcripts', subtitle = 'log-scale') +
  theme(axis.title = element_blank(), legend.position = 'none') +
  coord_flip()

# 필터링 전 발현된 유전자 수 (nFeature) 시각화 (선형 스케일)
p3 <- ggplot(cells, aes(x = sample_group, y = nFeature, fill = sample_group)) +
  geom_violin(draw_quantiles = c(0.5), scale = 'area', trim = FALSE) +
  geom_hline(yintercept = median_nFeature, color = 'black') +
  geom_hline(yintercept = thresholds_nFeature[2], color = 'red') +
  theme_bw() +
  scale_fill_manual(values = custom_colors$discrete) +
  scale_x_discrete(limits = rev(levels(cells$sample_group))) +
  scale_y_continuous(labels = scales::comma) +
  labs(title = 'Number of expressed genes', subtitle = 'linear scale') +
  theme(axis.title = element_blank(), legend.position = 'none') +
  coord_flip()

# 필터링 전 발현된 유전자 수 (nFeature) 시각화 (로그 스케일)
p4 <- ggplot(cells, aes(x = sample_group, y = nFeature, fill = sample_group)) +
  geom_violin(draw_quantiles = c(0.5), scale = 'area', trim = FALSE) +
  geom_hline(yintercept = median_nFeature, color = 'black') +
  geom_hline(yintercept = thresholds_nFeature[2], color = 'red') +
  theme_bw() +
  scale_fill_manual(values = custom_colors$discrete) +
  scale_x_discrete(limits = rev(levels(cells$sample_group))) +
  scale_y_log10(labels = scales::comma) +
  labs(title = 'Number of expressed genes', subtitle = 'log-scale') +
  theme(axis.title = element_blank(), legend.position = 'none') +
  coord_flip()

# 필터링 전 미토콘드리아 전사체 비율 (percent_MT) 시각화 (선형 스케일)
p5 <- ggplot(cells, aes(x = sample_group, y = percent_MT, fill = sample_group)) +
  geom_violin(draw_quantiles = c(0.5), scale = 'area', trim = FALSE) +
  geom_hline(yintercept = median_percent_MT, color = 'black') +
  geom_hline(yintercept = thresholds_percent_MT[2], color = 'red') +
  theme_bw() +
  scale_fill_manual(values = custom_colors$discrete) +
  scale_x_discrete(limits = rev(levels(cells$sample_group))) +
  scale_y_continuous(labels = scales::percent) +
  labs(title = 'Percent MT transcripts', subtitle = 'linear scale') +
  theme(axis.title = element_blank(), legend.position = 'none') +
  coord_flip()

# 필터링 임곗값 시각화 결과 저장
ggsave(
  'plots/qc_histogram_ncount_nfeature_percentmt_thresholds_PA_AC.png',
  p1 + p3 + p5 + p2 + p4 + plot_layout(ncol = 3),
  height = 7, width = 10
)


# nCount 대비 nFeature 분포 (임곗값 포함)
# PA와 AC 그룹에 맞게 nCount와 nFeature를 시각화
p <- ggplot(cells, aes(x = nCount, y = nFeature, color = percent_MT, shape = sample_group)) +
  geom_point(size = 0.5) +
  geom_hline(yintercept = thresholds_nFeature[2], color = 'red') + # nFeature 임곗값
  geom_vline(xintercept = thresholds_nCount[2], color = 'red') + # nCount 임곗값
  scale_x_continuous(name = 'Number of transcripts', labels = scales::comma) +
  scale_y_continuous(name = 'Number of expressed genes', labels = scales::comma) +
  theme_bw() +
  scale_color_viridis(
    name = 'Percent MT\ntranscripts',
    limits = c(0, 1),
    labels = scales::percent,
    guide = guide_colorbar(frame.colour = 'black', ticks.colour = 'black')
  ) +
  scale_shape_manual(values = c(16, 17)) + # PA와 AC 그룹에 다른 모양의 포인트 적용
  labs(title = 'nFeature vs nCount with thresholds for PA and AC') +
  theme(
    axis.title = element_blank(),
    legend.position = 'right',
    panel.grid.major.x = element_blank(),
    panel.grid.major.y = element_blank()
  )

# 필터링 기준이 적용된 nFeature vs nCount 분포 그래프 저장
ggsave('plots/qc_nfeature_over_ncount_thresholds_PA_AC.png', p, height = 4, width = 6)

# 3.1.5 필터링 후 데이터 확인 (After filtering)
# 필터링 기준을 적용하여 세포 선택 (PA와 AC 그룹으로 나누기)
cells_filtered <- cells %>%
  dplyr::filter(
    nCount >= thresholds_nCount[1], 
    nCount <= thresholds_nCount[2], 
    nFeature >= thresholds_nFeature[1], 
    nFeature <= thresholds_nFeature[2], 
    percent_MT >= thresholds_percent_MT[1], 
    percent_MT <= thresholds_percent_MT[2]
  )

# 필터링 후 5395개의 세포가 유지됨
# 각 샘플별 전사체 수, 발현된 유전자 수, 미토콘드리아 전사체 비율을 PA/AC 그룹으로 나누어 시각화

# 샘플별 그룹화 (PA와 AC로 나누기)
cells_filtered$sample_group <- ifelse(grepl("PA", cells_filtered$sample), "PA", "AC")

# 선형 스케일 전사체 수 (nCount)
p1 <- ggplot(cells_filtered, aes(x = sample_group, y = nCount, fill = sample_group)) +
  geom_violin(draw_quantiles = c(0.5), scale = 'area', trim = FALSE) +
  geom_hline(yintercept = median_nCount, color = 'black') +
  geom_hline(yintercept = thresholds_nCount[2], color = 'red') +
  theme_bw() +
  scale_fill_manual(values = custom_colors$discrete) +
  scale_x_discrete(limits = rev(levels(cells_filtered$sample_group))) +
  scale_y_continuous(labels = scales::comma) +
  labs(title = 'Number of transcripts', subtitle = 'linear scale') +
  theme(axis.title = element_blank(), legend.position = 'none') +
  coord_flip()

# 로그 스케일 전사체 수 (nCount)
p2 <- ggplot(cells_filtered, aes(x = sample_group, y = nCount, fill = sample_group)) +
  geom_violin(draw_quantiles = c(0.5), scale = 'area', trim = FALSE) +
  geom_hline(yintercept = median_nCount, color = 'black') +
  geom_hline(yintercept = thresholds_nCount[2], color = 'red') +
  theme_bw() +
  scale_fill_manual(values = custom_colors$discrete) +
  scale_x_discrete(limits = rev(levels(cells_filtered$sample_group))) +
  scale_y_log10(labels = scales::comma) +
  labs(title = 'Number of transcripts', subtitle = 'log-scale') +
  theme(axis.title = element_blank(), legend.position = 'none') +
  coord_flip()

# 선형 스케일 발현된 유전자 수 (nFeature)
p3 <- ggplot(cells_filtered, aes(x = sample_group, y = nFeature, fill = sample_group)) +
  geom_violin(draw_quantiles = c(0.5), scale = 'area', trim = FALSE) +
  geom_hline(yintercept = median_nFeature, color = 'black') +
  geom_hline(yintercept = thresholds_nFeature[2], color = 'red') +
  theme_bw() +
  scale_fill_manual(values = custom_colors$discrete) +
  scale_x_discrete(limits = rev(levels(cells_filtered$sample_group))) +
  scale_y_continuous(labels = scales::comma) +
  labs(title = 'Number of expressed genes', subtitle = 'linear scale') +
  theme(axis.title = element_blank(), legend.position = 'none') +
  coord_flip()

# 로그 스케일 발현된 유전자 수 (nFeature)
p4 <- ggplot(cells_filtered, aes(x = sample_group, y = nFeature, fill = sample_group)) +
  geom_violin(draw_quantiles = c(0.5), scale = 'area', trim = FALSE) +
  geom_hline(yintercept = median_nFeature, color = 'black') +
  geom_hline(yintercept = thresholds_nFeature[2], color = 'red') +
  theme_bw() +
  scale_fill_manual(values = custom_colors$discrete) +
  scale_x_discrete(limits = rev(levels(cells_filtered$sample_group))) +
  scale_y_log10(labels = scales::comma) +
  labs(title = 'Number of expressed genes', subtitle = 'log-scale') +
  theme(axis.title = element_blank(), legend.position = 'none') +
  coord_flip()

# 선형 스케일 미토콘드리아 전사체 비율 (percent_MT)
p5 <- ggplot(cells_filtered, aes(x = sample_group, y = percent_MT, fill = sample_group)) +
  geom_violin(draw_quantiles = c(0.5), scale = 'area', trim = FALSE) +
  geom_hline(yintercept = median_percent_MT, color = 'black') +
  geom_hline(yintercept = thresholds_percent_MT[2], color = 'red') +
  theme_bw() +
  scale_fill_manual(values = custom_colors$discrete) +
  scale_x_discrete(limits = rev(levels(cells_filtered$sample_group))) +
  scale_y_continuous(labels = scales::percent) +
  labs(title = 'Percent MT transcripts', subtitle = 'linear scale') +
  theme(axis.title = element_blank(), legend.position = 'none') +
  coord_flip()

# 필터링 후 nCount, nFeature, percent_MT의 분포 그래프 저장
ggsave(
  'plots/qc_histogram_ncount_nfeature_percentmt_filtered_PA_AC.png',
  p1 + p3 + p5 + p2 + p4 + plot_layout(ncol = 3),
  height = 7, width = 10
)

# 필터링 후 nFeature vs nCount 분포 (PA와 AC 그룹, 임곗값 포함)
p <- ggplot(cells_filtered, aes(nCount, nFeature, color = percent_MT, shape = sample_group)) +
  geom_point(size = 0.5) +
  geom_hline(yintercept = thresholds_nFeature[2], color = 'red') +  # nFeature 임곗값
  geom_vline(xintercept = thresholds_nCount[2], color = 'red') +  # nCount 임곗값
  scale_x_continuous(name = 'Number of transcripts', labels = scales::comma) +
  scale_y_continuous(name = 'Number of expressed genes', labels = scales::comma) +
  theme_bw() +
  scale_color_viridis(
    name = 'Percent MT\ntranscripts',
    limits = c(0, 1),
    labels = scales::percent,
    guide = guide_colorbar(frame.colour = 'black', ticks.colour = 'black')
  ) +
  labs(title = 'nFeature vs nCount with Thresholds') +  # 제목 추가
  theme(
    axis.title = element_blank(),
    legend.position = 'right'
  )

# 필터링 후 nFeature vs nCount 분포 그래프 저장
ggsave('plots/qc_nfeature_over_ncount_filtered_PA_AC.png', p, height = 4, width = 6)


# 필터링 후 유지할 세포 목록 저장
cells_to_keep <- cells_filtered$cell
length(cells_to_keep) # 필터링 후 유지된 세포 수: 10053개

# 3.2 유전자 필터링 (Filter genes)

# 각 유전자의 총 카운트 및 발현된 세포 수 계산
genes <- tibble(
  gene = rownames(transcripts$raw$merged), # 유전자 이름
  count = Matrix::rowSums(transcripts$raw$merged), # 총 카운트
  cells = Matrix::rowSums(transcripts$raw$merged != 0) # 발현된 세포 수
)

# 5개 이상의 세포에서 발현된 유전자만 선택
genes_to_keep <- genes %>% dplyr::filter(cells >= 5) %>% pull(gene)
length(genes_to_keep) # 16406개 유전자 유지

# 3.3 필터링된 전사체 행렬 생성 (Generate filtered transcript matrix)
transcripts$raw$filtered <- transcripts$raw$merged[genes_to_keep, cells_to_keep]

# 샘플별로 필터링 전후 세포 수 계산 (패턴 수정)
cells_per_sample_after_filtering <- tibble(
  sample = character(),
  before = numeric(),
  after = numeric()
)

for (i in c('AC', 'PA')) {
  tmp <- tibble(
    sample = i,
    before = grep(colnames(transcripts$raw$merged), pattern = paste0("Patient [1-3] ", i), value = FALSE) %>% length(),
    after = grep(colnames(transcripts$raw$filtered), pattern = paste0("Patient [1-3] ", i), value = FALSE) %>% length()
  )
  cells_per_sample_after_filtering <- bind_rows(cells_per_sample_after_filtering, tmp)
}

# 결과 확인
knitr::kable(cells_per_sample_after_filtering)

# ---------------------------
# 1) Seurat 객체 생성 (AC vs PA 비교를 위한 세팅)
seurat <- CreateSeuratObject(
  counts = transcripts$raw$filtered,  # 필터링된 데이터 사용
  min.cells = 0,
  min.features = 0
)

# 2) Seurat 메타데이터에서 샘플 정보를 'AC' 또는 'PA'로 설정
seurat@meta.data$sample <- Cells(seurat) %>%
  str_extract("Patient \\d+ (AC|PA)") %>%  # Patient 정보 제거 후 AC/PA만 남김
  str_extract("(AC|PA)")  # "Patient X " 제거 후 AC/PA 추출

# 3) sample을 factor로 변환하여 AC, PA 순서 유지
seurat@meta.data$sample <- factor(seurat@meta.data$sample, levels = c("AC", "PA"))

# 4) NA 값 제거 (샘플 정보가 없는 세포 제거)
seurat <- subset(seurat, cells = rownames(seurat@meta.data[!is.na(seurat@meta.data$sample), ]))

# 5) 샘플별 세포 수 확인 (AC vs PA)
temp_labels <- seurat@meta.data %>%
  group_by(sample) %>%
  tally()

# 6) 샘플별 nCount_RNA 시각화
p1 <- ggplot(seurat@meta.data, aes(x = sample, y = nCount_RNA, fill = sample)) +
  geom_violin(draw_quantiles = c(0.5), scale = 'area', trim = FALSE) +
  geom_half_boxplot(aes(x = sample, y = nCount_RNA, fill = sample),
                    side = 'r', outlier.color = NA, width = 0.4, show.legend = FALSE
  ) +
  geom_text(data = temp_labels, aes(x = sample, y = -Inf, label = paste0('n = ', format(n, big.mark = ',', trim = TRUE)), vjust = -1),
            color = 'black', size = 2.8
  ) +
  scale_fill_manual(values = custom_colors$discrete) +
  scale_y_continuous(labels = scales::comma, expand = c(0.08,0)) +
  theme_bw() +
  labs(x = '', y = 'Number of transcripts') +
  theme(panel.grid.major.x = element_blank(), axis.title.x = element_blank())

# 7) 샘플별 nFeature_RNA 시각화
p2 <- ggplot(seurat@meta.data, aes(x = sample, y = nFeature_RNA, fill = sample)) +
  geom_violin(draw_quantiles = c(0.5), scale = 'area', trim = FALSE) +
  geom_half_boxplot(aes(x = sample, y = nFeature_RNA, fill = sample),
                    side = 'r', outlier.color = NA, width = 0.4, show.legend = FALSE
  ) +
  geom_text(data = temp_labels, aes(x = sample, y = -Inf, label = paste0('n = ', format(n, big.mark = ',', trim = TRUE)), vjust = -1),
            color = 'black', size = 2.8
  ) +
  scale_fill_manual(values = custom_colors$discrete) +
  scale_y_continuous(name = 'Number of expressed genes', labels = scales::comma, expand = c(0.08,0)) +
  theme_bw() +
  theme(panel.grid.major.x = element_blank(), axis.title.x = element_blank())

# 시각화 결과 저장
ggsave(
  'plots/ncount_nfeature_by_sample_AC_PA.png',
  p1 + p2 + plot_layout(ncol = 2), height = 4, width = 10
)

# ----------------------------------

# 5 발현 데이터 정규화 (Normalize expression data)
seurat <- NormalizeData(
  seurat,
  assay = 'RNA',
  normalization.method = 'LogNormalize',
  scale.factor = median(seurat@meta.data$nCount_RNA)
)

# ----------------------------------

# 6 (선택) 데이터 통합 (Optional: Data integration using SCT + PCA)
seurat_list <- SplitObject(seurat, split.by = 'sample')

seurat_list <- lapply(
  X = seurat_list,
  FUN = function(x) {
    x <- SCTransform(x)
  }
)

seurat_features <- SelectIntegrationFeatures(
  seurat_list,
  nfeatures = 3000
)

seurat_list <- PrepSCTIntegration(
  seurat_list,
  anchor.features = seurat_features
)

seurat_list <- lapply(
  X = seurat_list,
  FUN = RunPCA,
  verbose = FALSE,
  features = seurat_features
)

seurat_anchors <- FindIntegrationAnchors(
  object.list = seurat_list,
  anchor.features = seurat_features,
  normalization.method = 'SCT',
  reduction = 'rpca'
)

seurat <- IntegrateData(
  anchorset = seurat_anchors,
  normalization.method = 'SCT'
)

# ----------------------------------

# 샘플 정보를 'AC' 또는 'PA'로 변환
seurat@meta.data$sample <- Cells(seurat) %>%
  str_extract("Patient \\d+ (AC|PA)") %>%  # 'Patient X ' 제거하고 AC/PA만 남김
  str_extract("(AC|PA)")  # AC 또는 PA 추출

# 샘플 정보를 factor로 변환하여 순서 유지
seurat@meta.data$sample <- factor(seurat@meta.data$sample, levels = c("AC", "PA"))

# NA 제거 (샘플 정보가 없는 세포 제거)
seurat <- subset(seurat, cells = rownames(seurat@meta.data[!is.na(seurat@meta.data$sample), ]))

# ---------------------------------------------

# 7. 주성분 분석 (Principal Component Analysis, PCA)
# 변수 유전자 설정 (FindVariableFeatures)
seurat <- FindVariableFeatures(seurat, selection.method = "vst", nfeatures = 2000)

# 변수 유전자 목록 확인
head(VariableFeatures(seurat))

# 주성분 분석 (Principal Component Analysis, PCA)
seurat <- RunPCA(seurat, assay = 'SCT', npcs = 50, features = VariableFeatures(seurat))

# 최적의 PCA 차원 결정
intrinsicDimension::maxLikGlobalDimEst(seurat@reductions$pca@cell.embeddings, k = 10)

# PCA 결과 시각화
p <- tibble(
  PC = 1:50,
  stdev = seurat@reductions$pca@stdev
) %>%
  ggplot(aes(PC, stdev)) +
  geom_point() +
  geom_vline(xintercept = 10, color = 'blue') +
  geom_vline(xintercept = 15, color = 'red') +
  theme_bw() +
  labs(x = 'Principal components', y = 'Standard deviation')

# PCA 결과 저장
ggsave('plots/principal_components.png', p, height = 4, width = 5)

# -----------------------------------------------------

# 8. 클러스터링 (Clustering)
seurat <- FindNeighbors(seurat, reduction = 'pca', dims = 1:15, graph.name = "pca_neighbors")
seurat <- FindClusters(seurat, resolution = 0.8, graph.name = "pca_neighbors")

# 클러스터별 세포 수 집계
temp_labels <- seurat@meta.data %>%
  group_by(seurat_clusters) %>%
  tally()

# 클러스터별 전사체 수(nCount_RNA) 시각화
p1 <- ggplot(seurat@meta.data, aes(seurat_clusters, nCount_RNA, fill = factor(seurat_clusters))) +
  geom_violin(trim = FALSE) +
  geom_boxplot(width = 0.2, outlier.shape = NA) +
  scale_fill_manual(values = custom_colors$discrete) +
  labs(y = 'Number of transcripts', x = 'Cluster') +
  theme_bw()

# 클러스터별 발현된 유전자 수(nFeature_RNA) 시각화
p2 <- ggplot(seurat@meta.data, aes(seurat_clusters, nFeature_RNA, fill = factor(seurat_clusters))) +
  geom_violin(trim = FALSE) +
  geom_boxplot(width = 0.2, outlier.shape = NA) +
  scale_fill_manual(values = custom_colors$discrete) +
  labs(y = 'Number of expressed genes', x = 'Cluster') +
  theme_bw()

# 결과 저장
ggsave('plots/ncount_nfeature_by_cluster.png', p1 + p2 + plot_layout(ncol = 1), height = 7, width = 14)

# -----------------------------------------------------

# 8.3 클러스터 안정성 평가 (Cluster stability)
library(SingleCellExperiment)
library(scater)
library(bluster)

# Seurat 객체에서 'SCT' assay의 'data' 레이어를 추출하여 변환
sce <- SingleCellExperiment(
  assays = list(logcounts = GetAssayData(seurat, assay = "SCT", layer = "data")),
  colData = seurat@meta.data
)

# Clustering 정보 추가
sce$seurat_clusters <- seurat$seurat_clusters

# PCA 실행
sce <- runPCA(sce, exprs_values = "logcounts", ncomponents = 50)

# PCA의 첫 15개 차원 사용
reducedDim(sce, 'PCA_sub') <- reducedDim(sce, 'PCA')[, 1:15, drop = FALSE]

# 부트스트랩 클러스터링 안정성 평가
ass_prob <- bootstrapStability(
  sce, 
  FUN = function(x) {
    g <- buildSNNGraph(x, use.dimred = 'PCA_sub')
    igraph::cluster_walktrap(g)$membership
  },
  clusters = sce$seurat_clusters
)

# 안정성 행렬 시각화
p <- as.data.frame(ass_prob) %>%
  rownames_to_column(var = 'cluster_1') %>%
  pivot_longer(cols = -cluster_1, names_to = 'cluster_2', values_to = 'probability') %>%
  ggplot(aes(cluster_2, cluster_1, fill = probability)) +
  geom_tile(color = 'white') +
  geom_text(aes(label = round(probability, digits = 2)), size = 2.5) +
  scale_fill_gradient(name = 'Probability', low = 'white', high = 'red', limits = c(0,1)) +
  theme_minimal()

# 결과 저장
ggsave('plots/cluster_stability.png', p, height = 6, width = 7)

# 8.4 실루엣 플롯 (Silhouette plot)

library(cluster)

# PCA의 첫 15개 차원을 사용하여 거리 행렬 계산
distance_matrix <- dist(Embeddings(seurat[['pca']])[, 1:15])

# 클러스터 레이블 추출
clusters <- seurat@meta.data$seurat_clusters

# 실루엣 점수 계산
silhouette <- silhouette(as.numeric(clusters), dist = distance_matrix)
seurat@meta.data$silhouette_score <- silhouette[,3] # 실루엣 점수 저장

# 평균 실루엣 점수 계산
mean_silhouette_score <- mean(seurat@meta.data$silhouette_score)

# 실루엣 점수 시각화
p <- seurat@meta.data %>%
  mutate(barcode = rownames(.)) %>%
  arrange(seurat_clusters, -silhouette_score) %>%
  mutate(barcode = factor(barcode, levels = barcode)) %>%
  ggplot() +
  geom_col(aes(barcode, silhouette_score, fill = seurat_clusters), show.legend = FALSE) +
  geom_hline(yintercept = mean_silhouette_score, color = 'red', linetype = 'dashed') + # 평균 실루엣 점수 표시
  scale_x_discrete(name = 'Cells') +
  scale_y_continuous(name = 'Silhouette score') +
  scale_fill_manual(values = custom_colors$discrete) +
  theme_bw() +
  theme(
    axis.title.x = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  )

# 실루엣 플롯 결과 저장
ggsave('plots/silhouette_plot.png', p, height = 4, width = 8)

# 8.5 클러스터 유사성 평가 (Cluster similarity)

library(bluster)
library(igraph)

# PCA의 첫 15개 차원을 사용하여 거리 행렬 계산
distance_matrix <- dist(reducedDim(sce, 'PCA_sub'))

# 최근접 이웃 그래프(SNN Graph) 생성
g <- buildSNNGraph(sce, use.dimred = 'PCA_sub')

# 클러스터 간 모듈러리티 비율 계산
modularity_ratio <- pairwiseModularity(
  g,  # 생성된 그래프 객체
  clusters = seurat$seurat_clusters  # 클러스터 정보
)

# 결과 확인
print(modularity_ratio)

# 로그 변환하여 시각화
ratio_to_plot <- log10(modularity_ratio + 1)

p <- ratio_to_plot %>%
  as_tibble() %>%
  rownames_to_column(var = 'cluster_1') %>%
  pivot_longer(
    cols = 2:ncol(.),
    names_to = 'cluster_2',
    values_to = 'probability'
  ) %>%
  mutate(
    cluster_1 = as.character(as.numeric(cluster_1) - 1),
    cluster_1 = factor(cluster_1, levels = rev(unique(cluster_1))),
    cluster_2 = factor(cluster_2, levels = unique(cluster_2))
  ) %>%
  ggplot(aes(cluster_2, cluster_1, fill = probability)) +
  geom_tile(color = 'white') +
  geom_text(aes(label = round(probability, digits = 2)), size = 2.5) +
  scale_x_discrete(name = 'Cluster', position = 'top') +
  scale_y_discrete(name = 'Cluster') +
  scale_fill_gradient(
    name = 'log10(ratio)', low = 'white', high = '#c0392b', na.value = '#bdc3c7',
    guide = guide_colorbar(
      frame.colour = 'black', ticks.colour = 'black', title.position = 'left',
      title.theme = element_text(hjust = 1, angle = 90),
      barwidth = 0.75, barheight = 10
    )
  ) +
  coord_fixed() +
  theme_bw() +
  theme(
    legend.position = 'right',
    panel.grid.major = element_blank()
  )

# 클러스터 유사성 결과 저장
ggsave('plots/cluster_similarity.png', p, height = 6, width = 7)

# 8.6 클러스터 트리 (Cluster tree)
seurat <- BuildClusterTree(
  seurat,
  dims = 1:15,
  reorder = FALSE,
  reorder.numeric = FALSE
)

# 트리 객체 가져오기
tree <- seurat@tools$BuildClusterTree
tree$tip.label <- paste0("Cluster ", tree$tip.label)

# 클러스터 계층적 트리 시각화
p <- ggtree::ggtree(tree, aes(x, y)) +
  scale_y_reverse() +
  ggtree::geom_tree() +
  ggtree::theme_tree() +
  ggtree::geom_tiplab(offset = 1) +
  ggtree::geom_tippoint(color = custom_colors$discrete[1:length(tree$tip.label)], shape = 16, size = 5) +
  coord_cartesian(clip = 'off') +
  theme(plot.margin = unit(c(0,2.5,0,0), 'cm'))

# 클러스터 트리 결과 저장
ggsave('plots/cluster_tree.png', p, height = 4, width = 6)

# 8.7 샘플과 클러스터 조성 비교 (Composition of samples and clusters)

# 각 샘플별 클러스터 분포 계산
table_samples_by_clusters <- seurat@meta.data %>%
  group_by(sample, seurat_clusters) %>%
  summarize(count = n(), .groups = 'drop') %>%
  spread(seurat_clusters, count, fill = 0) %>%
  ungroup() %>%
  mutate(total_cell_count = rowSums(select(., -sample))) %>%
  dplyr::select(sample, total_cell_count, everything()) %>%
  arrange(factor(sample, levels = levels(seurat@meta.data$sample)))

# 샘플별 클러스터 수 요약 테이블 출력
knitr::kable(table_samples_by_clusters)

# 각 클러스터별 샘플 분포 계산
table_clusters_by_samples <- seurat@meta.data %>%
  dplyr::rename('cluster' = 'seurat_clusters') %>%
  group_by(cluster, sample) %>%
  summarize(count = n(), .groups = 'drop') %>%
  spread(sample, count, fill = 0) %>%
  ungroup() %>%
  mutate(total_cell_count = rowSums(select(., -cluster))) %>%
  dplyr::select(cluster, total_cell_count, everything()) %>%
  arrange(factor(cluster, levels = levels(seurat@meta.data$seurat_clusters)))

# 클러스터별 샘플 수 요약 테이블 출력
knitr::kable(table_clusters_by_samples)

# 세포 수 분석 (Number of cells)

# 샘플별 총 세포 수 계산
temp_labels <- seurat@meta.data %>%
  group_by(sample) %>%
  tally()

# 샘플별 클러스터 구성 시각화 (세포 수 기준)
p1 <- table_samples_by_clusters %>%
  select(-c('total_cell_count')) %>%
  reshape2::melt(id.vars = 'sample') %>%
  mutate(sample = factor(sample, levels = levels(seurat@meta.data$sample))) %>%
  ggplot(aes(sample, value)) +
  geom_bar(aes(fill = variable), position = 'stack', stat = 'identity') +
  geom_text(
    data = temp_labels,
    aes(x = sample, y = Inf, label = paste0('n = ', format(n, big.mark = ',', trim = TRUE)), vjust = -1),
    color = 'black', size = 2.8
  ) +
  scale_fill_manual(name = 'Cluster', values = custom_colors$discrete) +
  scale_y_continuous(name = 'Number of cells', labels = scales::comma, expand = c(0.01, 0)) +
  coord_cartesian(clip = 'off') +
  theme_bw() +
  theme(
    legend.position = 'left',
    plot.title = element_text(hjust = 0.5),
    text = element_text(size = 16),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.title.x = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
    plot.margin = margin(t = 20, r = 0, b = 0, l = 0, unit = 'pt')
  )

# 클러스터별 총 세포 수 계산
temp_labels <- seurat@meta.data %>%
  group_by(seurat_clusters) %>%
  tally() %>%
  dplyr::rename('cluster' = seurat_clusters)

# 클러스터별 샘플 구성 시각화 (세포 수 기준)
p2 <- table_clusters_by_samples %>%
  select(-c('total_cell_count')) %>%
  reshape2::melt(id.vars = 'cluster') %>%
  mutate(cluster = factor(cluster, levels = levels(seurat@meta.data$seurat_clusters))) %>%
  ggplot(aes(cluster, value)) +
  geom_bar(aes(fill = variable), position = 'stack', stat = 'identity') +
  geom_text(
    data = temp_labels,
    aes(x = cluster, y = Inf, label = paste0('n = ', format(n, big.mark = ',', trim = TRUE)), vjust = -1),
    color = 'black', size = 2.8
  ) +
  scale_fill_manual(name = 'Sample', values = custom_colors$discrete) +
  scale_y_continuous(labels = scales::comma, expand = c(0.01, 0)) +
  coord_cartesian(clip = 'off') +
  theme_bw() +
  theme(
    legend.position = 'right',
    plot.title = element_text(hjust = 0.5),
    text = element_text(size = 16),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.title = element_blank(),
    plot.margin = margin(t = 20, r = 0, b = 0, l = 10, unit = 'pt')
  )

# 샘플별 및 클러스터별 세포 수 분포 시각화 결과 저장
ggsave(
  'plots/composition_samples_clusters_by_number.png',
  p1 + p2 +
    plot_layout(ncol = 2, widths = c(
      seurat@meta.data$sample %>% unique() %>% length(),
      seurat@meta.data$seurat_clusters %>% unique() %>% length()
    )),
  width = 18, height = 8
)

# 세포 비율 분석 (Percent of cells)

# 샘플별 총 세포 수 계산
temp_labels <- seurat@meta.data %>%
  group_by(sample) %>%
  tally()

# 샘플별 클러스터 구성 시각화 (비율 기준)
p1 <- table_samples_by_clusters %>%
  select(-c('total_cell_count')) %>%
  reshape2::melt(id.vars = 'sample') %>%
  mutate(sample = factor(sample, levels = levels(seurat@meta.data$sample))) %>%
  ggplot(aes(sample, value)) +
  geom_bar(aes(fill = variable), position = 'fill', stat = 'identity') +
  geom_text(
    data = temp_labels,
    aes(x = sample, y = Inf, label = paste0('n = ', format(n, big.mark = ',', trim = TRUE)), vjust = -1),
    color = 'black', size = 2.8
  ) +
  scale_fill_manual(name = 'Cluster', values = custom_colors$discrete) +
  scale_y_continuous(name = 'Percentage [%]', labels = scales::percent_format(), expand = c(0.01,0)) +
  coord_cartesian(clip = 'off') +
  theme_bw() +
  theme(
    legend.position = 'left',
    plot.title = element_text(hjust = 0.5),
    text = element_text(size = 16),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.title.x = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
    plot.margin = margin(t = 20, r = 0, b = 0, l = 0, unit = 'pt')
  )

# 클러스터별 총 세포 수 계산
temp_labels <- seurat@meta.data %>%
  group_by(seurat_clusters) %>%
  tally() %>%
  dplyr::rename('cluster' = seurat_clusters)

# 클러스터별 샘플 구성 시각화 (비율 기준)
p2 <- table_clusters_by_samples %>%
  select(-c('total_cell_count')) %>%
  reshape2::melt(id.vars = 'cluster') %>%
  mutate(cluster = factor(cluster, levels = levels(seurat@meta.data$seurat_clusters))) %>%
  ggplot(aes(cluster, value)) +
  geom_bar(aes(fill = variable), position = 'fill', stat = 'identity') +
  geom_text(
    data = temp_labels, aes(x = cluster, y = Inf, label = paste0('n = ', format(n, big.mark = ',', trim = TRUE)), vjust = -1),
    color = 'black', size = 2.8
  ) +
  scale_fill_manual(name = 'Sample', values = custom_colors$discrete) +
  scale_y_continuous(name = 'Percentage [%]', labels = scales::percent_format(), expand = c(0.01,0)) +
  coord_cartesian(clip = 'off') +
  theme_bw() +
  theme(
    legend.position = 'right',
    plot.title = element_text(hjust = 0.5),
    text = element_text(size = 16),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.title = element_blank(),
    plot.margin = margin(t = 20, r = 0, b = 0, l = 10, unit = 'pt')
  )

# 샘플별 및 클러스터별 세포 비율 시각화 결과 저장
ggsave(
  'plots/composition_samples_clusters_by_percent.png',
  p1 + p2 +
    plot_layout(ncol = 2, widths = c(
      seurat@meta.data$sample %>% unique() %>% length(),
      seurat@meta.data$seurat_clusters %>% unique() %>% length()
    )),
  width = 18, height = 8
)

# ----------------------------------

# 9. 세포 주기 분석 (Cell Cycle Analysis)

# 9.1 세포 주기 상태 추론 (Infer cell cycle status)
seurat <- CellCycleScoring(
  seurat,
  assay = 'SCT',
  s.features = cc.genes$s.genes,
  g2m.features = cc.genes$g2m.genes
)

# Seurat의 기본 Phase 컬럼을 사용자 정의로 변경
seurat@meta.data$cell_cycle_seurat <- seurat@meta.data$Phase
seurat@meta.data$Phase <- NULL
seurat@meta.data$cell_cycle_seurat <- factor(
  seurat@meta.data$cell_cycle_seurat, levels = c('G1', 'S', 'G2M')
)

# 9.2 샘플 및 클러스터별 세포 주기 구성 (Composition of samples and clusters by cell cycle)

# 샘플별 세포 주기 상태 분포
table_samples_by_cell_cycle <- seurat@meta.data %>%
  group_by(sample, cell_cycle_seurat) %>%
  summarize(count = n(), .groups = 'drop') %>%
  spread(cell_cycle_seurat, count, fill = 0) %>%
  ungroup() %>%
  mutate(total_cell_count = rowSums(select(., -sample))) %>%
  dplyr::select(sample, total_cell_count, everything()) %>%
  arrange(factor(sample, levels = levels(seurat@meta.data$cell_cycle_seurat)))

# 샘플별 세포 주기 상태 요약 테이블 출력
knitr::kable(table_samples_by_cell_cycle)

# 클러스터별 세포 주기 상태 분포
table_clusters_by_cell_cycle <- seurat@meta.data %>%
  dplyr::rename(cluster = seurat_clusters) %>%
  group_by(cluster, cell_cycle_seurat) %>%
  summarize(count = n(), .groups = 'drop') %>%
  spread(cell_cycle_seurat, count, fill = 0) %>%
  ungroup() %>%
  mutate(total_cell_count = rowSums(select(., -cluster))) %>%
  dplyr::select(cluster, total_cell_count, everything()) %>%
  arrange(factor(cluster, levels = levels(seurat@meta.data$cell_cycle_seurat)))

# 클러스터별 세포 주기 상태 요약 테이블 출력
knitr::kable(table_clusters_by_cell_cycle)

# 9.3 세포 주기 분포 시각화 (Number of cells by cell cycle)

# 샘플별 총 세포 수 계산
temp_labels <- seurat@meta.data %>%
  group_by(sample) %>%
  tally()

# 샘플별 세포 주기 상태 시각화 (막대 그래프, 세포 수 기준)
p1 <- table_samples_by_cell_cycle %>%
  select(-c('total_cell_count')) %>%
  reshape2::melt(id.vars = 'sample') %>%
  mutate(sample = factor(sample, levels = levels(seurat@meta.data$sample))) %>%
  ggplot(aes(sample, value)) +
  geom_bar(aes(fill = variable), position = 'stack', stat = 'identity') +
  geom_text(
    data = temp_labels,
    aes(x = sample, y = Inf, label = paste0('n = ', format(n, big.mark = ',', trim = TRUE)), vjust = -1),
    color = 'black', size = 2.8
  ) +
  scale_fill_manual(name = 'Cell cycle', values = custom_colors$cell_cycle) +
  scale_y_continuous(name = 'Number of cells', labels = scales::comma, expand = c(0.01,0)) +
  coord_cartesian(clip = 'off') +
  theme_bw() +
  theme(
    legend.position = 'none',
    plot.title = element_text(hjust = 0.5),
    text = element_text(size = 16),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.title.x = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
    plot.margin = margin(t = 20, r = 0, b = 0, l = 0, unit = 'pt')
  )

# 클러스터별 총 세포 수 계산
temp_labels <- seurat@meta.data %>%
  group_by(seurat_clusters) %>%
  tally() %>%
  dplyr::rename('cluster' = seurat_clusters)

# 클러스터별 세포 주기 상태 시각화 (막대 그래프, 세포 수 기준)
p2 <- table_clusters_by_cell_cycle %>%
  select(-c('total_cell_count')) %>%
  reshape2::melt(id.vars = 'cluster') %>%
  mutate(cluster = factor(cluster, levels = levels(seurat@meta.data$seurat_clusters))) %>%
  ggplot(aes(cluster, value)) +
  geom_bar(aes(fill = variable), position = 'stack', stat = 'identity') +
  geom_text(
    data = temp_labels,
    aes(x = cluster, y = Inf, label = paste0('n = ', format(n, big.mark = ',', trim = TRUE)), vjust = -1),
    color = 'black', size = 2.8
  ) +
  scale_fill_manual(name = 'Cell cycle', values = custom_colors$cell_cycle) +
  scale_y_continuous(name = 'Number of cells', labels = scales::comma, expand = c(0.01,0)) +
  coord_cartesian(clip = 'off') +
  theme_bw() +
  theme(
    legend.position = 'right',
    plot.title = element_text(hjust = 0.5),
    text = element_text(size = 16),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.title = element_blank(),
    plot.margin = margin(t = 20, r = 0, b = 0, l = 10, unit = 'pt')
  )

# 세포 주기 상태 시각화 결과 저장
ggsave(
  'plots/composition_samples_clusters_by_cell_cycle_by_number.png',
  p1 + p2 +
    plot_layout(ncol = 2, widths = c(
      seurat@meta.data$sample %>% unique() %>% length(),
      seurat@meta.data$seurat_clusters %>% unique() %>% length()
    )),
  width = 18, height = 8
)

# 9.4 세포 주기 비율 분석 (Percent of cells by cell cycle)

# 샘플별 총 세포 수 계산
temp_labels <- seurat@meta.data %>%
  group_by(sample) %>%
  tally()

# 샘플별 세포 주기 비율 시각화 (누적 막대 그래프, 비율 기준)
p1 <- table_samples_by_cell_cycle %>%
  select(-c('total_cell_count')) %>%
  reshape2::melt(id.vars = 'sample') %>%
  mutate(sample = factor(sample, levels = levels(seurat@meta.data$sample))) %>%
  ggplot(aes(sample, value)) +
  geom_bar(aes(fill = variable), position = 'fill', stat = 'identity') +
  geom_text(
    data = temp_labels,
    aes(x = sample, y = Inf, label = paste0('n = ', format(n, big.mark = ',', trim = TRUE)), vjust = -1),
    color = 'black', size = 2.8
  ) +
  scale_fill_manual(name = 'Cell cycle', values = custom_colors$cell_cycle) +
  scale_y_continuous(name = 'Percentage [%]', labels = scales::percent_format(), expand = c(0.01,0)) +
  coord_cartesian(clip = 'off') +
  theme_bw() +
  theme(
    legend.position = 'none',
    plot.title = element_text(hjust = 0.5),
    text = element_text(size = 16),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.title.x = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
    plot.margin = margin(t = 20, r = 0, b = 0, l = 0, unit = 'pt')
  )

# 클러스터별 총 세포 수 계산
temp_labels <- seurat@meta.data %>%
  group_by(seurat_clusters) %>%
  tally() %>%
  dplyr::rename('cluster' = seurat_clusters)

# 클러스터별 세포 주기 비율 시각화 (누적 막대 그래프, 비율 기준)
p2 <- table_clusters_by_cell_cycle %>%
  select(-c('total_cell_count')) %>%
  reshape2::melt(id.vars = 'cluster') %>%
  mutate(cluster = factor(cluster, levels = levels(seurat@meta.data$seurat_clusters))) %>%
  ggplot(aes(cluster, value)) +
  geom_bar(aes(fill = variable), position = 'fill', stat = 'identity') +
  geom_text(
    data = temp_labels,
    aes(x = cluster, y = Inf, label = paste0('n = ', format(n, big.mark = ',', trim = TRUE)), vjust = -1),
    color = 'black', size = 2.8
  ) +
  scale_fill_manual(name = 'Cell cycle', values = custom_colors$cell_cycle) +
  scale_y_continuous(name = 'Percentage [%]', labels = scales::percent_format(), expand = c(0.01,0)) +
  coord_cartesian(clip = 'off') +
  theme_bw() +
  theme(
    legend.position = 'right',
    plot.title = element_text(hjust = 0.5),
    text = element_text(size = 16),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.title = element_blank(),
    plot.margin = margin(t = 20, r = 0, b = 0, l = 10, unit = 'pt')
  )

# 세포 주기 상태 비율 시각화 결과 저장
ggsave(
  'plots/composition_samples_clusters_by_cell_cycle_by_percent.png',
  p1 + p2 +
    plot_layout(ncol = 2, widths = c(
      seurat@meta.data$sample %>% unique() %>% length(),
      seurat@meta.data$seurat_clusters %>% unique() %>% length()
    )),
  width = 18, height = 8
)

# ------------------------------------------

# 10. 차원 축소 분석 (Dimensional Reduction)

# 10.1 SNN 그래프 (Shared Nearest Neighbor Graph)
library(ggnetwork)

seurat <- FindNeighbors(seurat, reduction = "pca", dims = 1:30)


# Seurat 객체의 SNN 그래프 변환
SCT_snn <- seurat@graphs$SCT_snn %>%
  as.matrix() %>%
  ggnetwork() %>%
  left_join(seurat@meta.data %>% mutate(vertex.names = rownames(.)), by = 'vertex.names')

# 샘플별 SNN 그래프 시각화
p1 <- ggplot(SCT_snn, aes(x = x, y = y, xend = xend, yend = yend)) +
  geom_edges(color = 'grey50', alpha = 0.05) +
  geom_nodes(aes(color = sample), size = 0.5) +
  scale_color_manual(
    name = 'Sample', values = custom_colors$discrete,
    guide = guide_legend(ncol = 1, override.aes = list(size = 2))
  ) +
  theme_blank() +
  theme(legend.position = 'left') +
  annotate(
    geom = 'text', x = Inf, y = -Inf,
    label = paste0('n = ', format(nrow(seurat@meta.data), big.mark = ',', trim = TRUE)),
    vjust = -1.5, hjust = 1.25, color = 'black', size = 2.5
  )

# 클러스터별 SNN 그래프 시각화
p2 <- ggplot(SCT_snn, aes(x = x, y = y, xend = xend, yend = yend)) +
  geom_edges(color = 'grey50', alpha = 0.05) +
  geom_nodes(aes(color = seurat_clusters), size = 0.5) +
  scale_color_manual(
    name = 'Cluster', values = custom_colors$discrete,
    guide = guide_legend(ncol = 1, override.aes = list(size = 2))
  ) +
  theme_blank() +
  theme(legend.position = 'right') +
  annotate(
    geom = 'text', x = Inf, y = -Inf,
    label = paste0('n = ', format(nrow(seurat@meta.data), big.mark = ',', trim = TRUE)),
    vjust = -1.5, hjust = 1.25, color = 'black', size = 2.5
  )

# SNN 그래프 저장
ggsave(
  'plots/snn_graph_by_sample_cluster.png',
  p1 + p2 + plot_layout(ncol = 2),
  height = 5, width = 11
)

# 10.2 UMAP (Uniform Manifold Approximation and Projection)

# 10.2.1 UMAP 계산
seurat <- RunUMAP(
  seurat,
  reduction.name = 'UMAP',
  reduction = 'pca',
  dims = 1:15,
  n.components = 2,
  seed.use = 100
)

# 10.2.2 UMAP 개요

# UMAP 상에서 nCount_RNA 분포
plot_umap_by_nCount <- bind_cols(seurat@meta.data, as.data.frame(seurat@reductions$UMAP@cell.embeddings)) %>%
  ggplot(aes(UMAP_1, UMAP_2, color = nCount_RNA)) +
  geom_point(size = 0.2) +
  theme_bw() +
  scale_color_viridis(
    guide = guide_colorbar(frame.colour = 'black', ticks.colour = 'black'),
    labels = scales::comma
  ) +
  labs(color = 'Number of\ntranscripts') +
  theme(legend.position = 'left') +
  coord_fixed()

# UMAP 상에서 샘플 분포
plot_umap_by_sample <- bind_cols(seurat@meta.data, as.data.frame(seurat@reductions$UMAP@cell.embeddings)) %>%
  ggplot(aes(UMAP_1, UMAP_2, color = sample)) +
  geom_point(size = 0.2) +
  theme_bw() +
  scale_color_manual(values = custom_colors$discrete) +
  labs(color = 'Sample') +
  guides(colour = guide_legend(override.aes = list(size = 2))) +
  theme(legend.position = 'right') +
  coord_fixed()

# UMAP 상에서 클러스터 분포
plot_umap_by_cluster <- bind_cols(seurat@meta.data, as.data.frame(seurat@reductions$UMAP@cell.embeddings)) %>%
  ggplot(aes(UMAP_1, UMAP_2, color = seurat_clusters)) +
  geom_point(size = 0.2) +
  theme_bw() +
  scale_color_manual(
    name = 'Cluster', values = custom_colors$discrete,
    guide = guide_legend(ncol = 2, override.aes = list(size = 2))
  ) +
  theme(legend.position = 'left') +
  coord_fixed()

# UMAP 상에서 세포 주기 상태 분포
plot_umap_by_cell_cycle <- bind_cols(seurat@meta.data, as.data.frame(seurat@reductions$UMAP@cell.embeddings)) %>%
  ggplot(aes(UMAP_1, UMAP_2, color = cell_cycle_seurat)) +
  geom_point(size = 0.2) +
  theme_bw() +
  scale_color_manual(values = custom_colors$cell_cycle) +
  labs(color = 'Cell cycle') +
  guides(colour = guide_legend(override.aes = list(size = 2))) +
  theme(legend.position = 'right') +
  coord_fixed()

# UMAP 개요 그래프 저장
ggsave(
  'plots/umap.png',
  plot_umap_by_nCount + plot_umap_by_sample +
    plot_umap_by_cluster + plot_umap_by_cell_cycle +
    plot_layout(ncol = 2),
  height = 6,
  width = 8.5
)

# 10.2.3 샘플별 UMAP 분포

# 샘플별 총 세포 수 계산
temp_labels <- seurat@meta.data %>%
  group_by(sample) %>%
  tally()

# 샘플별 UMAP 시각화
p <- bind_cols(seurat@meta.data, as.data.frame(seurat@reductions$UMAP@cell.embeddings)) %>%
  ggplot(aes(UMAP_1, UMAP_2, color = sample)) +
  geom_point(size = 0.2, show.legend = FALSE) +
  geom_text(
    data = temp_labels,
    aes(x = Inf, y = -Inf, label = paste0('n = ', format(n, big.mark = ',', trim = TRUE)), vjust = -1.5, hjust = 1.25),
    color = 'black', size = 2.8
  ) +
  theme_bw() +
  scale_color_manual(values = custom_colors$discrete) +
  coord_fixed() +
  theme(
    legend.position = 'none',
    strip.text = element_text(face = 'bold', size = 10)
  ) +
  facet_wrap(~sample, ncol = 5)

# 샘플별 UMAP 분포 저장
ggsave(
  'plots/umap_by_sample_split.png',
  p,
  height = 4,
  width = 10
)

# 10.2.4 샘플별(PA, AC) 클러스터 중심점 계산 및 UMAP 시각화
library(ggplot2)
library(dplyr)

# 샘플별 데이터 분리 (AC, PA)
seurat_list <- SplitObject(seurat, split.by = "sample")

# UMAP 컬럼명 확인
print(head(colnames(seurat_list$AC@reductions$UMAP@cell.embeddings)))

# 샘플별 UMAP 시각화 및 클러스터 중심점 계산 함수 정의
plot_umap_per_sample <- function(seurat_obj, sample_name) {
  
  # UMAP 좌표 컬럼명을 자동으로 가져오기
  umap_data <- as.data.frame(seurat_obj@reductions$UMAP@cell.embeddings)
  colnames(umap_data) <- c("UMAP_1", "UMAP_2")  # UMAP 차원 이름을 통일
  
  # 클러스터 중심 좌표 계산
  UMAP_centers <- tibble(
    UMAP_1 = umap_data$UMAP_1,
    UMAP_2 = umap_data$UMAP_2,
    cluster = seurat_obj@meta.data$seurat_clusters
  ) %>%
    group_by(cluster) %>%
    summarize(x = median(UMAP_1), y = median(UMAP_2))
  
  # UMAP 시각화
  p <- bind_cols(seurat_obj@meta.data, umap_data) %>%
    ggplot(aes(UMAP_1, UMAP_2, color = seurat_clusters)) +
    geom_point(size = 0.2) +
    geom_label(
      data = UMAP_centers,
      mapping = aes(x, y, label = cluster),
      size = 4.5,
      fill = 'white',
      color = 'black',
      fontface = 'bold',
      alpha = 0.5,
      label.size = 0,
      show.legend = FALSE
    ) +
    theme_bw() +
    scale_color_manual(
      name = 'Cluster', values = custom_colors$discrete,
      guide = guide_legend(override.aes = list(size = 2))
    ) +
    theme(legend.position = 'right') +
    coord_fixed() +
    annotate(
      geom = 'text', x = Inf, y = -Inf,
      label = paste0('n = ', format(nrow(seurat_obj@meta.data), big.mark = ',', trim = TRUE)),
      vjust = -1.5, hjust = 1.25, color = 'black', size = 2.5
    ) +
    labs(title = paste("UMAP by Clusters (", sample_name, ")", sep = ""))
  
  # 결과 저장
  ggsave(
    filename = paste0('plots/umap_by_clusters_', sample_name, '.png'),
    plot = p, height = 5.5, width = 5.5
  )
}

# AC와 PA 각각에 대해 UMAP 플롯 생성
plot_umap_per_sample(seurat_list$AC, "AC")
plot_umap_per_sample(seurat_list$PA, "PA")


# 10.3 Poincaré 맵 (현재 작업 중)
# 자세한 내용은 블로그 "Animating the dynamics of training a Poincaré map" 참고

# -----------------------------------

# 11 SingleR을 활용한 세포 유형 주석 (Cell Type Annotation)
library(SingleR)
library(Seurat)
library(ggplot2)
library(dplyr)

# 11.1 세포 유형 할당 (Assign labels)

# BlueprintEncodeData()를 이용한 참조 데이터 로드
singler_ref <- BlueprintEncodeData()

# Seurat 객체를 AC와 PA로 분리
seurat_list <- SplitObject(seurat, split.by = "sample")

# 각 샘플별로 SingleR 분석을 수행하는 함수 정의
run_singleR_per_sample <- function(seurat_obj, sample_name) {
  
  # 무작위성 고정
  set.seed(32)  # 어떤 숫자든 괜찮지만, 재현을 위해 고정
  
  # SingleR을 사용하여 세포 유형 예측 수행
  singler_results <- SingleR::SingleR(
    test = GetAssayData(seurat_obj, assay = 'SCT', slot = 'data'),
    ref = singler_ref,
    labels = singler_ref@colData@listData$label.main
  )
  
  # 결과 저장
  saveRDS(singler_results, paste0("data/singler_results_", sample_name, ".rds"))
  
  # SingleR 결과 히트맵 생성
  p <- plotScoreHeatmap(
    singler_results,
    show.labels = TRUE,
    annotation_col = data.frame(
      donor = seurat_obj@meta.data$sample,
      row.names = rownames(singler_results)
    )
  )
  
  # 히트맵 저장
  ggsave(
    paste0('plots/singler_score_heatmap_', sample_name, '.png'), p,
    height = 11, width = 14
  )
  
  # Seurat 객체의 메타데이터에 예측된 세포 유형 추가
  seurat_obj@meta.data$cell_type_singler <- singler_results@listData$labels
  
  # 세포 유형 예측 점수 저장
  singler_scores <- singler_results@listData$scores %>%
    as_tibble() %>%
    dplyr::mutate(assigned_score = NA)
  
  for (i in seq_len(nrow(singler_scores))) {
    singler_scores$assigned_score[i] <- singler_scores[[singler_results@listData$labels[i]]][i]
  }
  
  seurat_obj@meta.data$cell_type_singler_score <- singler_scores$assigned_score
  
  return(seurat_obj)
}

# AC와 PA 각각에 대해 SingleR 수행
seurat_list$AC <- run_singleR_per_sample(seurat_list$AC, "AC")
seurat_list$PA <- run_singleR_per_sample(seurat_list$PA, "PA")





#####
library(SingleR)
library(SingleCellExperiment)

# BlueprintEncodeData 참조 데이터 로드
singler_ref <- BlueprintEncodeData()

# Seurat 객체를 SingleCellExperiment로 변환
seurat_sce <- as.SingleCellExperiment(seurat, assay = "SCT")

# SingleCellExperiment에서 'data' 슬롯 가져오기
sce_data <- assay(seurat_sce, "logcounts")  # 'data' 슬롯은 logcounts에 해당

# SingleR 실행
set.seed(32)
singler_results <- SingleR(
  test = sce_data,
  ref = singler_ref,
  labels = singler_ref$label.main
)

# 결과를 Seurat 객체의 meta.data에 추가
seurat$cell_type_singler <- singler_results$labels
seurat$cell_type_singler_score <- apply(singler_results$scores, 1, max)  # 최고 점수를 저장

# 결과 확인
table(seurat$cell_type_singler)










# 11.2 샘플별(PA, AC) 클러스터 및 세포 유형별 할당 점수 (Assignment Score)
library(ggplot2)
library(dplyr)
library(scales)
library(viridis)

# Seurat 객체를 AC와 PA로 분리
seurat_list <- SplitObject(seurat, split.by = "sample")

# 샘플별 할당 점수 시각화 함수 정의
plot_assignment_scores <- function(seurat_obj, sample_name) {
  
  # NA 값 제거 및 데이터 필터링
  seurat_obj@meta.data <- seurat_obj@meta.data %>%
    filter(!is.na(cell_type_singler_score))  # NA 제거
  
  # 샘플별 세포 수 계산
  temp_labels <- seurat_obj@meta.data %>%
    group_by(sample, .drop = FALSE) %>%
    tally()
  
  # 샘플별 할당 점수 분포 시각화
  p1 <- ggplot(seurat_obj@meta.data, aes(x = sample, y = cell_type_singler_score, fill = sample)) +
    geom_violin(trim = FALSE, show.legend = FALSE) +
    geom_boxplot(width = 0.2, outlier.shape = NA, show.legend = FALSE) +
    geom_text(data = temp_labels, aes(x = sample, y = -Inf, label = paste0('n = ', format(n, big.mark = ',', trim = TRUE)), vjust = -1), color = 'black', size = 2.8) +
    scale_fill_manual(values = c("PA" = "blue", "AC" = "red")) +
    scale_y_continuous(name = 'Assignment score', labels = scales::comma, limits = c(0,1)) +
    theme_bw() +
    labs(title = paste("Samples (", sample_name, ")", sep = "")) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))
  
  # 클러스터별 세포 수 계산
  temp_labels <- seurat_obj@meta.data %>%
    group_by(seurat_clusters, .drop = FALSE) %>%
    tally()
  
  # 클러스터별 할당 점수 분포 시각화 (색상 자동 확장)
  p2 <- ggplot(seurat_obj@meta.data, aes(x = factor(seurat_clusters), y = cell_type_singler_score, fill = factor(seurat_clusters))) +
    geom_violin(trim = FALSE, show.legend = FALSE) +
    geom_boxplot(width = 0.2, outlier.shape = NA, show.legend = FALSE) +
    geom_text(data = temp_labels, aes(x = factor(seurat_clusters), y = -Inf, label = paste0('n = ', format(n, big.mark = ',', trim = TRUE)), vjust = -1), color = 'black', size = 2.8) +
    scale_fill_viridis(discrete = TRUE, option = "D") +  # 자동 색상 확장
    scale_y_continuous(name = 'Assignment score', labels = scales::comma, limits = c(0,1)) +
    labs(title = paste("Clusters (", sample_name, ")", sep = "")) +
    theme_bw()
  
  # 세포 유형별 세포 수 계산
  temp_labels <- seurat_obj@meta.data %>%
    group_by(cell_type_singler, .drop = FALSE) %>%
    tally()
  
  # 세포 유형별 할당 점수 분포 시각화 (색상 자동 확장)
  p3 <- ggplot(seurat_obj@meta.data, aes(x = cell_type_singler, y = cell_type_singler_score, fill = cell_type_singler)) +
    geom_violin(trim = FALSE, show.legend = FALSE) +
    geom_boxplot(width = 0.2, outlier.shape = NA, show.legend = FALSE) +
    geom_text(data = temp_labels, aes(x = cell_type_singler, y = -Inf, label = paste0('n = ', format(n, big.mark = ',', trim = TRUE)), vjust = -1), color = 'black', size = 2.8) +
    scale_fill_viridis(discrete = TRUE, option = "C") +  # 자동 색상 확장
    scale_y_continuous(name = 'Assignment score', labels = scales::comma, limits = c(0,1)) +
    labs(title = paste("Cell types (", sample_name, ")", sep = "")) +
    theme_bw()
  
  # 결과 저장
  ggsave(
    filename = paste0('plots/singler_assignment_scores_', sample_name, '.png'),
    plot = p1 + p2 + p3 + plot_layout(ncol = 1), height = 13, width = 16
  )
}

# 데이터 필터링: 최소 3개 이상 존재하는 세포 유형만 유지
seurat_list$PA@meta.data <- seurat_list$PA@meta.data %>%
  group_by(cell_type_singler, .drop = FALSE) %>%
  filter(n() > 2) %>%
  ungroup()

seurat_list$AC@meta.data <- seurat_list$AC@meta.data %>%
  group_by(cell_type_singler, .drop = FALSE) %>%
  filter(n() > 2) %>%
  ungroup()

# AC와 PA 각각에 대해 할당 점수 시각화
plot_assignment_scores(seurat_list$AC, "AC")
plot_assignment_scores(seurat_list$PA, "PA")








# 11.3 샘플 및 클러스터별 세포 유형 구성 분석
library(ggplot2)
library(dplyr)
library(tidyr)
library(scales)
library(reshape2)
library(patchwork)

# Seurat 객체 내 세포 유형 컬럼 확인
valid_cell_type_col <- "cell_type_singler"

if (!valid_cell_type_col %in% colnames(seurat@meta.data)) {
  stop("Error: SingleR 세포 주석 컬럼이 존재하지 않습니다. SingleR을 다시 실행하세요.")
}

# 샘플별 세포 유형 빈도 테이블 생성
table_samples_by_cell_type <- seurat@meta.data %>%
  dplyr::group_by(sample, .data[[valid_cell_type_col]]) %>%
  dplyr::summarize(count = n(), .groups = "drop") %>%
  tidyr::spread(.data[[valid_cell_type_col]], count, fill = 0) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(total_cell_count = rowSums(.[, -1])) %>%
  dplyr::select(sample, total_cell_count, everything())

# 샘플별 세포 유형 분포 출력
table_samples_by_cell_type %>% knitr::kable()

# 클러스터별 세포 유형 빈도 테이블 생성
table_clusters_by_cell_type <- seurat@meta.data %>%
  dplyr::group_by(seurat_clusters, .data[[valid_cell_type_col]]) %>%
  dplyr::summarize(count = n(), .groups = "drop") %>%
  tidyr::spread(.data[[valid_cell_type_col]], count, fill = 0) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(total_cell_count = rowSums(.[, -1])) %>%
  dplyr::select(seurat_clusters, total_cell_count, everything()) %>%
  dplyr::rename(cluster = seurat_clusters)

# 클러스터별 세포 유형 분포 출력
table_clusters_by_cell_type %>% knitr::kable()

# ----------------------------
# 샘플별 세포 유형 시각화
temp_labels <- seurat@meta.data %>%
  group_by(sample) %>%
  tally()

p1 <- table_samples_by_cell_type %>%
  select(-total_cell_count) %>%
  reshape2::melt(id.vars = "sample") %>%
  mutate(sample = factor(sample, levels = unique(seurat@meta.data$sample))) %>%
  ggplot(aes(sample, value)) +
  geom_bar(aes(fill = variable), position = "stack", stat = "identity", show.legend = FALSE) +
  geom_text(
    data = temp_labels,
    aes(
      x = sample,
      y = Inf,
      label = paste0("n = ", format(n, big.mark = ",", trim = TRUE)),
      vjust = -1
    ),
    color = "black", size = 2.8
  ) +
  scale_fill_viridis_d(name = "Cell type") +
  scale_y_continuous(
    name = "Number of cells",
    labels = scales::comma,
    expand = c(0.01, 0)
  ) +
  coord_cartesian(clip = "off") +
  theme_bw() +
  theme(
    legend.position = "none",
    plot.title = element_text(hjust = 0.5),
    text = element_text(size = 16),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.title.x = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
    plot.margin = margin(t = 20, r = 0, b = 0, l = 0, unit = "pt")
  )

# ----------------------------
# 클러스터별 세포 유형 시각화
temp_labels <- seurat@meta.data %>%
  group_by(seurat_clusters) %>%
  tally() %>%
  dplyr::rename(cluster = seurat_clusters)

p2 <- table_clusters_by_cell_type %>%
  select(-total_cell_count) %>%
  reshape2::melt(id.vars = "cluster") %>%
  mutate(cluster = factor(cluster, levels = unique(seurat@meta.data$seurat_clusters))) %>%
  ggplot(aes(cluster, value)) +
  geom_bar(aes(fill = variable), position = "stack", stat = "identity") +
  geom_text(
    data = temp_labels,
    aes(
      x = cluster,
      y = Inf,
      label = paste0("n = ", format(n, big.mark = ",", trim = TRUE)),
      vjust = -1
    ), color = "black", size = 2.8
  ) +
  scale_fill_viridis_d(name = "Cell type") +
  scale_y_continuous(
    name = "Number of cells",
    labels = scales::comma,
    expand = c(0.01, 0)
  ) +
  coord_cartesian(clip = "off") +
  theme_bw() +
  theme(
    legend.position = "right",
    plot.title = element_text(hjust = 0.5),
    text = element_text(size = 16),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.title = element_blank(),
    plot.margin = margin(t = 20, r = 0, b = 0, l = 10, unit = "pt")
  )




# 11.3.1 올루비얼(흐름) 플롯을 이용한 샘플과 클러스터 간의 관계 시각화

# 샘플 및 클러스터 이름 가져오기
samples <- levels(seurat@meta.data$sample)
clusters <- levels(seurat@meta.data$seurat_clusters)

# 샘플 및 클러스터별 색상 할당
color_assignments <- setNames(
  c(custom_colors$discrete[1:length(samples)], custom_colors$discrete[1:length(clusters)]),
  c(samples, clusters)
)

# 데이터 준비 (샘플과 클러스터 간의 세포 수)
data <- seurat@meta.data %>%
  group_by(sample, seurat_clusters) %>%
  tally() %>%
  ungroup() %>%
  gather_set_data(1:2) %>%
  dplyr::mutate(
    x = factor(x, levels = unique(x)),
    y = factor(y, levels = unique(y))
  )

# 샘플 및 클러스터 레이블 생성 (위치 조정)
data_labels <- tibble(
  group = c(
    rep('sample', length(samples)),
    rep('seurat_clusters', length(clusters))
  )
) %>%
  mutate(
    hjust = ifelse(group == 'sample', 1, 0),
    nudge_x = ifelse(group == 'sample', -0.1, 0.1)
  )

# 올루비얼(흐름) 플롯 생성
p1 <- ggplot(data, aes(x, id = id, split = y, value = n)) +
  geom_parallel_sets(aes(fill = seurat_clusters), alpha = 0.75, axis.width = 0.15) +
  geom_parallel_sets_axes(aes(fill = y), color = 'black', axis.width = 0.1) +
  geom_text(
    aes(y = n, split = y), stat = 'parallel_sets_axes', fontface = 'bold',
    hjust = data_labels$hjust, nudge_x = data_labels$nudge_x
  ) +
  scale_x_discrete(labels = c('Sample', 'Cluster')) +
  scale_fill_manual(values = color_assignments) +
  theme_bw() +
  theme(
    legend.position = 'none',
    axis.title = element_blank(),
    axis.text.x = element_text(face = 'bold', colour = 'black', size = 15),
    axis.text.y = element_blank(),
    axis.ticks = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank()
  )

# 그래프 저장
ggsave('plots/alluvial_plot_samples_clusters.png', p1, width = 8, height = 5)

# 클러스터별 세포 유형 분석 및 시각화

# 클러스터 및 세포 유형 목록 설정
clusters <- levels(seurat@meta.data$seurat_clusters)
cell_types <- sort(unique(seurat@meta.data$cell_type_singler))

# 색상 매핑
color_assignments <- setNames(
  c(custom_colors$discrete[1:length(clusters)], custom_colors$discrete[1:length(cell_types)]),
  c(clusters, cell_types)
)

# 데이터 정리 및 변환 (클러스터 -> 세포 유형 흐름)
data <- seurat@meta.data %>%
  dplyr::rename(cell_type = cell_type_singler) %>%
  dplyr::mutate(cell_type = factor(cell_type, levels = cell_types)) %>%
  group_by(seurat_clusters, cell_type) %>%
  tally() %>%
  ungroup() %>%
  gather_set_data(1:2) %>%
  dplyr::mutate(
    x = factor(x, levels = unique(x)),
    y = factor(y, levels = c(clusters, cell_types))
  )

# 축 레이블 설정
data_labels <- tibble(
  group = c(
    rep('seurat_clusters', length(clusters)),
    rep('cell_type', length(cell_types))
  )
) %>%
  mutate(
    hjust = ifelse(group == 'seurat_clusters', 1, 0),
    nudge_x = ifelse(group == 'seurat_clusters', -0.1, 0.1)
  )

# 클러스터 -> 세포 유형 올루비얼(흐름) 플롯 생성
p2 <- ggplot(data, aes(x, id = id, split = y, value = n)) +
  geom_parallel_sets(aes(fill = seurat_clusters), alpha = 0.75, axis.width = 0.15) +
  geom_parallel_sets_axes(aes(fill = y), color = 'black', axis.width = 0.1) +
  geom_text(
    aes(y = n, split = y), stat = 'parallel_sets_axes', fontface = 'bold',
    hjust = data_labels$hjust, nudge_x = data_labels$nudge_x
  ) +
  scale_x_discrete(labels = c('Cluster', 'Cell type')) +
  scale_fill_manual(values = color_assignments) +
  theme_bw() +
  theme(
    legend.position = 'none',
    axis.title = element_blank(),
    axis.text.x = element_text(face = 'bold', colour = 'black', size = 15),
    axis.text.y = element_blank(),
    axis.ticks = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank()
  )

# 샘플 <-> 클러스터 플롯과 클러스터 <-> 세포 유형 플롯 결합 및 저장
ggsave(
  'plots/samples_clusters_cell_types_alluvial.png',
  p1 + p2 + plot_layout(ncol = 2),
  height = 6, width = 8
)


# 필요한 패키지 로드
library(dplyr)
library(tibble)
library(circlize)

# 1. 메타데이터 tibble로 변환 및 클러스터 factor 처리
meta <- as_tibble(seurat@meta.data) %>%
  mutate(
    seurat_clusters = as.factor(seurat_clusters),
    cell_type_singler = as.factor(cell_type_singler)
  )

# 2. 샘플 ↔ 클러스터 관계 데이터
data_samples_clusters <- meta %>%
  group_by(sample, seurat_clusters) %>%
  tally(name = "n") %>%
  ungroup() %>%
  dplyr::rename(from = sample, to = seurat_clusters)

# 3. 클러스터 ↔ 세포유형 관계 데이터
data_clusters_celltypes <- meta %>%
  group_by(seurat_clusters, cell_type_singler) %>%
  tally(name = "n") %>%
  ungroup() %>%
  dplyr::rename(from = seurat_clusters, to = cell_type_singler)

# 4. 두 관계 데이터를 하나로 통합
chord_data <- bind_rows(data_samples_clusters, data_clusters_celltypes)

# 5. 모든 항목(샘플, 클러스터, 세포유형)에 색상 지정
groups <- unique(c(chord_data$from, chord_data$to))
color_palette <- setNames(rainbow(length(groups)), groups)

# 6. 서큘러 플롯 그리기
circos.clear()
chordDiagram(
  x = chord_data,
  grid.col = color_palette,
  transparency = 0.3,
  annotationTrack = "grid",
  preAllocateTracks = list(track.height = 0.05)
)

# 7. 섹터 레이블 정렬
circos.trackPlotRegion(
  track.index = 1,
  panel.fun = function(x, y) {
    sector_name <- get.cell.meta.data("sector.index")
    circos.text(
      x = mean(get.cell.meta.data("xlim")),
      y = 0,
      labels = sector_name,
      facing = "clockwise",
      niceFacing = TRUE,
      adj = c(0, 0.5),
      cex = 0.75
    )
  },
  bg.border = NA
)

# 8. 저장
png("plots/chord_plot_sample_cluster_celltype.png", width = 1000, height = 1000, res = 150)
chordDiagram(chord_data, grid.col = color_palette, transparency = 0.3)
dev.off()


# Seurat 객체 이름: seurat (가정)
# OASL 발현 기준으로 그룹 나누기 및 DEG 수행

### 1. AC / PA 정보 → Disease_Status로 변환 ------------------------------
seurat$Disease_Status <- ifelse(seurat$sample == "PA", "NonDisease", "Disease")


### 2. OASL 발현 추출 및 그룹 분리 -------------------------------------
oasl_expr <- FetchData(seurat, vars = "OASL")$OASL
names(oasl_expr) <- colnames(seurat)

# 셀 구분
ctl_cells <- colnames(seurat)[seurat$Disease_Status == "NonDisease"]
dis_cells <- colnames(seurat)[seurat$Disease_Status == "Disease"]

# 중앙값 기준으로 나눔
ctl_median <- median(oasl_expr[ctl_cells], na.rm = TRUE)
dis_median <- median(oasl_expr[dis_cells], na.rm = TRUE)

# 그룹 지정
seurat$OASL_group <- NA
seurat$OASL_group[ctl_cells] <- ifelse(oasl_expr[ctl_cells] >= ctl_median, "High", "Low")
seurat$OASL_group[dis_cells] <- ifelse(oasl_expr[dis_cells] >= dis_median, "High", "Low")


### 3. 네 개의 세포 그룹 정의 ------------------------------------------
cells_disease_high <- colnames(seurat)[seurat$Disease_Status == "Disease" & seurat$OASL_group == "High"]
cells_disease_low  <- colnames(seurat)[seurat$Disease_Status == "Disease" & seurat$OASL_group == "Low"]
cells_control_high <- colnames(seurat)[seurat$Disease_Status == "NonDisease" & seurat$OASL_group == "High"]
cells_control_low  <- colnames(seurat)[seurat$Disease_Status == "NonDisease" & seurat$OASL_group == "Low"]


### 4. DEG 분석 실행 ---------------------------------------------------
# 필요시 logfc.threshold, min.pct 조정 가능
deg1 <- FindMarkers(seurat, ident.1 = cells_disease_high, ident.2 = cells_disease_low,
                    logfc.threshold = 0, min.pct = 0.1)
deg2 <- FindMarkers(seurat, ident.1 = cells_control_high, ident.2 = cells_control_low,
                    logfc.threshold = 0, min.pct = 0.1)
deg3 <- FindMarkers(seurat, ident.1 = cells_disease_high, ident.2 = cells_control_high,
                    logfc.threshold = 0, min.pct = 0.1)
deg4 <- FindMarkers(seurat, ident.1 = cells_disease_low, ident.2 = cells_control_low,
                    logfc.threshold = 0, min.pct = 0.1)


### 5. DEG 결과 저장 ----------------------------------------------------
write.csv(deg1, "scRNA_DEG1_Disease_High_vs_Low.csv")
write.csv(deg2, "scRNA_DEG2_NonDisease_High_vs_Low.csv")
write.csv(deg3, "scRNA_DEG3_High_Disease_vs_NonDisease.csv")
write.csv(deg4, "scRNA_DEG4_Low_Disease_vs_NonDisease.csv")





# ░▒▓ 패키지 로드 ▓▒░
library(Seurat)
library(dplyr)
library(ggplot2)
library(patchwork)
library(RColorBrewer)

# ░▒▓ Step 1: UMAP + 메타데이터 병합 ▓▒░
umap_df <- bind_cols(
  seurat@meta.data,
  as.data.frame(Embeddings(seurat, "UMAP"))
)

# ░▒▓ Step 2: CellType 색상 정의 ▓▒░
celltype_levels <- unique(umap_df$cell_type_singler)
custom_colors <- list()
custom_colors$celltype <- setNames(
  RColorBrewer::brewer.pal(n = max(3, min(length(celltype_levels), 12)), name = "Set1") %>%
    rep(length.out = length(celltype_levels)),
  celltype_levels
)

# ░▒▓ Step 3: Cell Type Label 중심 계산 ▓▒░
label_df_non <- umap_df %>%
  filter(Disease_Status == "NonDisease") %>%
  group_by(cell_type_singler) %>%
  summarise(UMAP_1 = median(UMAP_1), UMAP_2 = median(UMAP_2), .groups = "drop")

label_df_dis <- umap_df %>%
  filter(Disease_Status == "Disease") %>%
  group_by(cell_type_singler) %>%
  summarise(UMAP_1 = median(UMAP_1), UMAP_2 = median(UMAP_2), .groups = "drop")

# ░▒▓ Step 4: A - CellType UMAP (NonDisease + label) ▓▒░
p1 <- umap_df %>%
  filter(Disease_Status == "NonDisease") %>%
  ggplot(aes(UMAP_1, UMAP_2, color = cell_type_singler)) +
  geom_point(size = 0.2, alpha = 0.9) +
  geom_text(data = label_df_non, aes(label = cell_type_singler),
            color = "black", size = 3.2, fontface = "bold", check_overlap = TRUE) +
  scale_color_manual(values = custom_colors$celltype) +
  coord_fixed() +
  theme_bw() +
  labs(title = "A. Cell Type Annotation (NonDisease)") +
  theme(
    legend.position = "right",
    plot.title = element_text(face = "bold")
  )

# ░▒▓ Step 5: B - CellType UMAP (Disease + label) ▓▒░
p2 <- umap_df %>%
  filter(Disease_Status == "Disease") %>%
  ggplot(aes(UMAP_1, UMAP_2, color = cell_type_singler)) +
  geom_point(size = 0.2, alpha = 0.9) +
  geom_text(data = label_df_dis, aes(label = cell_type_singler),
            color = "black", size = 3.2, fontface = "bold", check_overlap = TRUE) +
  scale_color_manual(values = custom_colors$celltype) +
  coord_fixed() +
  theme_bw() +
  labs(title = "B. Cell Type Annotation (Disease)") +
  theme(
    legend.position = "right",
    plot.title = element_text(face = "bold")
  )

# ░▒▓ Step 6: OASL 발현 시각화 함수 ▓▒░
oasl_color_map <- function(data) {
  data$OASL_for_plot <- ifelse(data$OASL > 0, data$OASL, NA)
  
  ggplot() +
    geom_point(data = data %>% filter(OASL <= 0),
               aes(x = UMAP_1, y = UMAP_2),
               color = "lightgrey", size = 0.25, alpha = 0.6) +
    geom_point(data = data %>% filter(OASL > 0),
               aes(x = UMAP_1, y = UMAP_2, color = OASL_for_plot),
               size = 0.25, alpha = 0.9) +
    scale_color_gradient(low = "pink", high = "red", name = "OASL") +
    coord_fixed() +
    theme_bw()
}

# ░▒▓ Step 7: C/D - OASL Expression UMAP (그대로 유지) ▓▒░
p3 <- umap_df %>%
  filter(Disease_Status == "NonDisease") %>%
  oasl_color_map() +
  labs(title = "C. OASL Expression (NonDisease)") +
  theme(
    plot.title = element_text(face = "bold"),
    legend.position = "right"
  )

p4 <- umap_df %>%
  filter(Disease_Status == "Disease") %>%
  oasl_color_map() +
  labs(title = "D. OASL Expression (Disease)") +
  theme(
    plot.title = element_text(face = "bold"),
    legend.position = "right"
  )

# ░▒▓ Step 8: 2x2 Figure ▓▒░
final_plot <- (p1 | p2) / (p3 | p4)

# ░▒▓ Step 9: 저장 ▓▒░
ggsave("plots/umap_2x2_Celltype_OASL_labeled.png", final_plot, width = 12, height = 9)






# ░▒▓ 패키지 로드 ▓▒░
library(Seurat)
library(dplyr)
library(ggplot2)
library(patchwork)
library(RColorBrewer)
library(tidyr)  # complete 사용을 위해 필요

# ░▒▓ Step 1: UMAP + 메타데이터 병합 ▓▒░
umap_df <- bind_cols(
  seurat@meta.data,
  as.data.frame(Embeddings(seurat, "UMAP"))
)

# ░▒▓ Step 2: Cell Type 순서 정의 (뒤집은 순서) ▓▒░
cell_order <- rev(c(
  "CD4+ T-cells", "CD8+ T-cells", "B-cells", "NK cells", "DC", "Monocytes", "Macrophages",
  "Endothelial cells", "Fibroblasts", "Smooth muscle",
  "Chondrocytes", "Adipocytes", "Astrocytes", "Skeletal muscle", "HSC", "Myocytes"
))

# ░▒▓ Step 3: Cell Type + Condition별 평균 OASL ▓▒░
avg_expr <- umap_df %>%
  filter(cell_type_singler %in% cell_order) %>%
  group_by(Disease_Status, cell_type_singler) %>%
  summarise(mean_OASL = mean(OASL, na.rm = TRUE), .groups = "drop")

# ░▒▓ Step 4: 유의성 검정 (t-test) ▓▒░
pvals <- umap_df %>%
  filter(cell_type_singler %in% cell_order) %>%
  group_by(cell_type_singler) %>%
  summarise(p = tryCatch(
    t.test(OASL ~ Disease_Status)$p.value,
    error = function(e) NA_real_
  ), .groups = "drop") %>%
  mutate(sig = case_when(
    is.na(p) ~ "",
    p < 0.001 ~ "***",
    p < 0.01 ~ "**",
    p < 0.05 ~ "*",
    TRUE ~ ""
  ))

# ░▒▓ Step 5: Z-score 계산 및 병합 ▓▒░
avg_expr <- avg_expr %>%
  group_by(Disease_Status) %>%
  mutate(OASL_z = scale(mean_OASL)[, 1]) %>%
  ungroup() %>%
  complete(cell_type_singler, Disease_Status, fill = list(mean_OASL = NA, OASL_z = NA)) %>%
  left_join(pvals, by = "cell_type_singler") %>%
  mutate(
    cell_type_singler = factor(cell_type_singler, levels = cell_order),
    Disease_Status = factor(Disease_Status, levels = c("NonDisease", "Disease")),
    label = paste0(sprintf("%.2f", OASL_z), sig)
  )

# ░▒▓ Step 6: Heatmap 시각화 ▓▒░
p_heatmap <- ggplot(avg_expr, aes(x = Disease_Status, y = cell_type_singler, fill = OASL_z)) +
  geom_tile(color = "black", linewidth = 0.3) +  # 테두리 색상과 두께 조절
  geom_text(aes(label = label), size = 3.2) +
  scale_fill_gradient2(low = "skyblue", mid = "white", high = "red", midpoint = 0,
                       name = "Z-score\n(OASL)") +
  theme_minimal() +
  labs(title = "OASL Expression (Z-score) by Cell Type & Condition",
       x = NULL, y = NULL) +
  theme(
    axis.text.x = element_text(face = "bold"),
    axis.text.y = element_text(face = "bold"),
    panel.grid = element_blank(),  # 그리드 제거
    plot.title = element_text(face = "bold", hjust = 0.5)
  )


# ░▒▓ Step 7: 저장 ▓▒░
ggsave("plots/heatmap_OASL_by_Celltype_DiseaseStatus_Zscore_pval.png", p_heatmap, width = 5.5, height = 7)


# 저장 코드: 분석 데이터 & Seurat 객체 저장
#  저장 디렉토리 생성 (없으면 자동 생성) 
dir.create("data", showWarnings = FALSE)

# 1. Seurat 객체 저장 (UMAP, 클러스터링 포함 전체 상태) 
saveRDS(seurat, file = "data/seurat_OASL_processed.rds")

# 2. UMAP + 메타데이터 테이블 저장 (Z-score 분석 포함) 
saveRDS(umap_df, file = "data/umap_df_OASL_analysis.rds")

#  3. 평균 발현 + Z-score + 유의성 테이블 저장 (heatmap용) 
saveRDS(avg_expr, file = "data/avg_expr_OASL_heatmap.rds")

# (선택) CSV 파일도 함께 저장 (외부 공유용) ▓▒░
# write.csv(umap_df, "data/umap_df_OASL_analysis.csv", row.names = FALSE)
# write.csv(avg_expr, "data/avg_expr_OASL_heatmap.csv", row.names = FALSE)


# 불러오기
seurat <- readRDS("data/seurat_OASL_processed.rds")
umap_df <- readRDS("data/umap_df_OASL_analysis.rds")
avg_expr <- readRDS("data/avg_expr_OASL_heatmap.rds")

# CoDEG1 에 대해 scRNA-seq에서 검증
# ░▒▓ 패키지 로드 ▓▒░
library(Seurat)
library(ggplot2)
library(patchwork)

# ░▒▓ 1. 허브 유전자 리스트 ▓▒░
hub_genes <- c("PTPRC", "CCL5", "LCP2", "TLR2", "CD3D", "CD8A", "LCK", "IKZF1", "TRAF3IP3", "PTPN6")

# ░▒▓ 2. 세포 유형 기준 UMAP 발현 시각화 ▓▒░
Idents(seurat) <- "cell_type_singler"

FeaturePlot(seurat, features = hub_genes, reduction = "UMAP", cols = c("gray90", "red"), pt.size = 0.3)

# ░▒▓ 3. DotPlot (세포 유형별 발현 및 비율) ▓▒░
DotPlot(seurat, features = hub_genes, group.by = "cell_type_singler") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# ░▒▓ 4. RidgePlot: OASL 그룹별 발현 분포 (전체 세포) ▓▒░
Idents(seurat) <- "OASL_group"

RidgePlot(seurat, features = hub_genes, group.by = "OASL_group", ncol = 2)

# ░▒▓ 5. Macrophage에서만 OASL High vs Low 비교 (VlnPlot) ▓▒░
seurat_mac <- subset(seurat, subset = cell_type_singler == "Macrophages")

VlnPlot(seurat_mac, features = hub_genes, group.by = "OASL_group", pt.size = 0.2) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


# CoDEG 3-4에 대해 검증
# ─────────────────────────────────────────────────────────────
# CoDEG3 & CoDEG4 허브 유전자 시각화 
# ─────────────────────────────────────────────────────────────

library(Seurat)
library(ggplot2)
library(patchwork)
library(dplyr)

# 1. Seurat 객체 identity 설정
Idents(seurat) <- "cell_type_singler"

# 2. 허브 유전자 리스트
hub_list <- list(
  CoDEG3 = c("OASL","SKAP1","TYROBP","CD14","CCL5","LIPA","FBLN1",
             "CD3G","HLA-DMA","CD5","RASAL3","THEMIS"),
  CoDEG4 = c("OASL","HLA-DMA","HLA-DRB1","CTSB","TYROBP","COL5A1",
             "PCOLCE","HLA-DMB","CCL3","CD74","IL7R","CD83")
)

# 3. 면역세포 타입 정의
immune_types <- c("CD4+ T-cells","CD8+ T-cells","NK cells",
                  "B-cells","DC","Monocytes","Macrophages")

# 4. OASL_group 컬럼 확인
if (!"OASL_group" %in% colnames(seurat@meta.data)) {
  stop("메타데이터에 'OASL_group' 컬럼이 없습니다.")
}

for (cond in names(hub_list)) {
  hub_genes <- hub_list[[cond]]
  
  # ──────────────────────────
  # A) UMAP: 여러 패널에 제목 붙이기
  umap_plots <- FeaturePlot(
    seurat,
    features  = hub_genes,
    reduction  = "UMAP",
    cols       = c("gray90","red"),
    pt.size    = 0.3,
    combine    = FALSE
  )
  p_umap <- wrap_plots(umap_plots, ncol = min(4, length(umap_plots))) +
    plot_annotation(title = paste0(cond, " UMAP")) &
    theme(plot.title = element_text(face="bold", hjust=0.5))
  ggsave(paste0(cond, "_UMAP.png"), p_umap, width=12, height=6, dpi=300)
  
  # ──────────────────────────
  # B) DotPlot
  p_dot <- DotPlot(
    seurat,
    features = hub_genes,
    group.by = "cell_type_singler"
  ) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    ggtitle(paste0(cond, " DotPlot")) +
    theme(plot.title = element_text(face="bold", hjust=0.5))
  ggsave(paste0(cond, "_DotPlot.png"), p_dot, width=8, height=6, dpi=300)
  
  # ──────────────────────────
  # C) RidgePlot: OASL High vs Low
  Idents(seurat) <- "OASL_group"
  ridge_plots <- RidgePlot(
    seurat,
    features = hub_genes,
    group.by = "OASL_group",
    ncol     = 2,
    combine  = FALSE
  )
  p_ridge <- wrap_plots(ridge_plots, ncol=2) +
    plot_annotation(title = paste0(cond, " RidgePlot (OASL High vs Low)")) &
    theme(plot.title = element_text(face="bold", hjust=0.5))
  ggsave(paste0(cond, "_RidgePlot.png"), p_ridge, width=12, height=8, dpi=300)
  
  # ──────────────────────────
  # D) 면역세포 ViolinPlot
  seu_immune <- subset(seurat, subset = cell_type_singler %in% immune_types)
  Idents(seu_immune) <- "OASL_group"
  violin_list <- VlnPlot(
    seu_immune,
    features  = hub_genes,
    pt.size   = 0,
    group.by  = "OASL_group",
    split.by  = "cell_type_singler",
    combine   = FALSE
  )
  p_violin <- wrap_plots(violin_list, ncol=2) +
    plot_annotation(title = paste0(cond, " Immune Cells Violin (OASL High vs Low)")) &
    theme(plot.title = element_text(face="bold", hjust=0.5))
  ggsave(paste0(cond, "_Immune_Violin.png"), p_violin, width=12, height=10, dpi=300)
  
  # ──────────────────────────
  # identity 복원
  Idents(seurat) <- "cell_type_singler"
}
