setwd("/path/to/work/ZS.Sigmoid/")

system("cp /work/ZS/CCLE.4VAE.ZS.tsv /work/ZS.Sigmoid/")
system("cp /work/ZS/TCGA.4VAE.ZS.tsv /work/ZS.Sigmoid/")
dir.create("/work/ZS.Sigmoid/result")
dir.create("/work/ZS.Sigmoid/result.EN/")
dir.create("/work/ZS.Sigmoid/result.EN/dr.CCLE")
dir.create("/work/ZS.Sigmoid/result.EN/dr.GDSC")

for(b in 1:100){
	cmd = paste("python3 run.VAE_100_0005_100_100.Sigmoid.py ", b, 
	" /work/ZS.Sigmoid/ CCLE.4VAE.ZS.tsv TCGA.4VAE.ZS.tsv CCLE.latent.tsv CCLE.weight.tsv TCGA.latent.tsv encoder.hdf5 decoder.hdf5", sep="" )
	system(cmd)
}
