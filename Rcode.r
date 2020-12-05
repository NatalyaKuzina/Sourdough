#загружаем пакеты и данные (https://benjjneb.github.io/dada2/dada-installation.html)
library(dada2); packageVersion("dada2")
library(vegan)
library("ggplot2"); packageVersion("ggplot2")
library("ape")
library(phyloseq); packageVersion("phyloseq")
library(ShortRead)
packageVersion("ShortRead")
library(Biostrings); packageVersion("Biostrings")
vignette("phyloseq-basics")

path <- "~/storage/VNIIGRJ/leaven/leaven-2/leaven_1_16S/"  #меняем на свой путь к данным (не забываем про экранирование слэша для windows  \\)
list.files(path)
#читаем имена файлов fastq и выполняем некоторые манипуляции со строками, чтобы получить совпадающие списки прямого и обратного файлов fastq.
fnFs <- sort(list.files(path, pattern="_R1_001.fastq", full.names = TRUE))
fnRs <- sort(list.files(path, pattern="_R2_001.fastq", full.names = TRUE))
#получаем имена образцов
sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 1)

#начнем с визуализации профилей качества форвард-ридов
plotQualityProfile(fnFs[6:7])
#визуализируем реверсы
plotQualityProfile(fnRs[6:7])

#назначим имена файлов для отфильтрованных файлов fastq.gz
filtFs <- file.path(path, "filtered", paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(path, "filtered", paste0(sample.names, "_R_filt.fastq.gz"))
names(filtFs) <- sample.names
names(filtRs) <- sample.names
#обрезаем и фильтруем, определяем truncLen исходя из визуализации качества ридов и некоторых личных предпочтений к обработке данных
out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, truncLen=c(240,180),
                     maxN=0, maxEE=c(2,2), truncQ=2, rm.phix=TRUE,
                     compress=TRUE, multithread=TRUE) # On Windows set multithread=FALSE
#проверка, что все отработало верно
head(out)
#Алгоритм DADA2 использует модель параметрической ошибки(err), и каждый набор данных ампликона имеет свой набор ошибок. 
#learnErrorsМетод узнает эту модель ошибки из данных, с помощью переменного оценивания частоты ошибок и 
#вывода образца композиции, пока они не сходятся на совместно последовательное решение. 
#Алгоритм начинается с первоначального предположения, для которого используются максимально возможные коэффициенты ошибок в этих данных 
#(коэффициенты ошибок, если верна только наиболее распространенная последовательность, а все остальные ошибки).
errF <- learnErrors(filtFs, multithread=TRUE)
errR <- learnErrors(filtRs, multithread=TRUE)
#в качестве проверки работоспособности визуализируем оценочные коэффициенты ошибок:
plotErrors(errF, nominalQ=TRUE)

#вывод образцов с высоким разрешением по данным об ампликонах Illumina
dadaFs <- dada(filtFs, err=errF, multithread=TRUE)
dadaRs <- dada(filtRs, err=errR, multithread=TRUE)
#просматриваем возвращенный dada-class объект:
dadaFs[[1]]

#объединяем прямое и обратное чтение вместе
mergers <- mergePairs(dadaFs, filtFs, dadaRs, filtRs, verbose=TRUE)
#проверяем наш дата-фрейм просмотрев первые образцы
head(mergers[[1]])

#строим таблицу вариантов последовательности ампликонов (ASV), версию таблицы OTU с более высоким разрешением, созданную традиционными методами
seqtab <- makeSequenceTable(mergers)
#проверяем размерность
dim(seqtab)

#проверяем распределение длинн последовательностей
table(nchar(getSequences(seqtab)))
head((nchar(getSequences(seqtab))))
#удаляем химеры
seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)
#еще пара проверок того, что все работает верно
dim(seqtab.nochim)
sum(seqtab.nochim)/sum(seqtab)

#смотрим конечные результаты числа прочтений в образцах, оценивая, не было ли сильных потерь на каком-либо из этапов:
getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN), rowSums(seqtab.nochim))
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
rownames(track) <- sample.names
head(track)
write.csv(track, "track.csv")

#назначаем таксономию вариантам последовательности 
#Пакет DADA2 предоставляет для этой цели встроенную реализацию метода наивного байесовского классификатора. 
#assignTaxonomyФункция принимает в качестве входных данных набор последовательностей, 
#чтобы быть классифицированы и обучающего набора эталонных последовательностей с известной систематики и 
#выводит таксономические назначения, по меньшей мере minBoot.
taxa <- assignTaxonomy(seqtab.nochim, "~/storage/tax/silva_nr_v132_train_set.fa.gz", multithread=TRUE)
taxa <- addSpecies(taxa, "~/storage/tax/silva_species_assignment_v132.fa.gz")
taxa.print <- taxa #данную таблицу можно посмотреть в окне справа (environment), нажав на иконку в виде таблички
rownames(taxa.print) <- NULL
#проверка
head(taxa.print)
#выводим готовую табличку в csv файл
write.csv(taxa.print, "taxa.print.csv") 
