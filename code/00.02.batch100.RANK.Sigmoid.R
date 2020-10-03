setwd("/path/to/work/RANK.Sigmoid/")

system("cp /work/RANK/CCLE.4VAE.RANK.tsv /work/RANK.Sigmoid/")
system("cp /work/RANK/TCGA.4VAE.RANK.tsv /work/RANK.Sigmoid/")
dir.create("/work/RANK.Sigmoid/result")
dir.create("/work/RANK.Sigmoid/result.EN/")
dir.create("/work/RANK.Sigmoid/result.EN/dr.CCLE")
dir.create("/work/RANK.Sigmoid/result.EN/dr.GDSC")

for(b in 1:100){
	cmd = paste("python3 run.VAE_100_0005_100_100.Sigmoid.py ", b, 
	" /work/RANK.Sigmoid/ CCLE.4VAE.RANK.tsv TCGA.4VAE.RANK.tsv CCLE.latent.tsv CCLE.weight.tsv TCGA.latent.tsv encoder.hdf5 decoder.hdf5", sep="" )
	system(cmd)
}
