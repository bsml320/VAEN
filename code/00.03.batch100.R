setwd("/work/")

### VAE models will be saved in the folder result
dir.create("result")

for(b in 1:100){
	cmd = paste("python3 run.VAE_100_0005_100_100.py ", b, " /work/ CCLE.4VAE-peer.tsv TCGA.4VAE-peer.tsv GSE65185.4VAE-peer.tsv CCLE.latent.tsv CCLE.weight.tsv TCGA.latent.tsv GSE65185.latent.tsv encoder.hdf5 decoder.hdf5", sep="" )
	system(cmd)
}
