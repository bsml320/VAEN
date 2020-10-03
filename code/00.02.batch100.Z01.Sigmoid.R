setwd("/path/to/work/Z01.Sigmoid/")

system("cp /work/Z01/CCLE.4VAE.Z01.tsv /work/Z01.Sigmoid/")
system("cp /work/Z01/TCGA.4VAE.Z01.tsv /work/Z01.Sigmoid/")
dir.create("/work/Z01.Sigmoid/result")
dir.create("/work/Z01.Sigmoid/result.EN/")
dir.create("/work/Z01.Sigmoid/result.EN/dr.CCLE")
dir.create("/work/Z01.Sigmoid/result.EN/dr.GDSC")

for(b in 1:100){
	cmd = paste("python3 run.VAE_100_0005_100_100.Sigmoid.py ", b, 
	" /work/Z01.Sigmoid/ CCLE.4VAE.Z01.tsv TCGA.4VAE.Z01.tsv CCLE.latent.tsv CCLE.weight.tsv TCGA.latent.tsv encoder.hdf5 decoder.hdf5", sep="" )
	system(cmd)
}
