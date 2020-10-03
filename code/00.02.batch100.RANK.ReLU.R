setwd("/path/to/work/RANK.ReLU/")

system("cp /work/RANK/CCLE.4VAE.RANK.tsv /work/RANK.ReLU/")
system("cp /work/RANK/TCGA.4VAE.RANK.tsv /work/RANK.ReLU/")
dir.create("/work/RANK.ReLU/result")
dir.create("/work/RANK.ReLU/result.EN/")
dir.create("/work/RANK.ReLU/result.EN/dr.CCLE")
dir.create("/work/RANK.ReLU/result.EN/dr.GDSC")

for(b in 1:100){
	cmd = paste("python3 run.VAE_100_0005_100_100.ReLU.py ", b, 
	" /work/RANK.ReLU/ CCLE.4VAE.RANK.tsv TCGA.4VAE.RANK.tsv CCLE.latent.tsv CCLE.weight.tsv TCGA.latent.tsv encoder.hdf5 decoder.hdf5", sep="" )
	system(cmd)
}
