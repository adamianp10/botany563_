Course project

My research aims to reconstruct a complete phylogeny of the neotropical clades of the orchid genus Vanilla (subgen. Vanilla and sect. Xanata). Since I don't have my own data yet, I'll use accessions posted in GenBank, the DNA regions chosen are matK and rbcL.

#git commands
  git add . 
  git commit –m “uploaded README.txt” 
  git push 

##ALIGNMENT-CLUSTALW

I had several problems trying to install clustalw. Finally I used this command and it worked!

	conda create -n clustalw2 -c biobuilds -y clustalw

To activate clustalw the command is:

	conda activate clustalw2

To run the alignment first you need to be located inside the folder that has the fasta file, then run the following command:

	clustalw2 -ALIGN -INFILE=primatesAA.fasta -OUTFILE=primatesAA-aligned.fasta -OUTPUT=PHYLIP

To merge various .fasta files use "cat" command. In this case "> 1.fasta" means that all sequences will me merge on that file, this will overwrite whatever is inside that file

	cat 2.fasta 3.fasta 4.fasta 5.fasta 6.fasta 7.fasta 8.fasta 9.fasta 10.fasta 11.fasta 12.fasta 13.fasta 14.fasta 15.fasta 16.fasta 17.fasta 18.fasta 19.fasta 20.fasta 21.fasta >1.fasta

[March 2]

I had problems installing R, this error pop up:

Error reading R script (), system error 2 (No such file or directory); Unable to find libR.dylib in expected locationswithin R Home directory /Library/Frameworks/R.framework/Resources.

I solved installing first the R-4.1.pkg from here: https://cran.r-project.org/bin/macosx/

Later I tried to install these packages

install.packages("adegenet", dep=TRUE)
install.packages("phangorn", dep=TRUE)

But this error keeps popping up

Warning in install.packages :
  installation of package ‘ape’ had non-zero exit status

I (apparently) solved it using tools/install packages/ape on R and responding with a NO to the following question:

Do you want to install from sources the package which needs compilation? (Yes/no/cancel) 

Now the packages run without any problems

library(ape)
library(adegenet)
library(phangorn)


#[March 4]

#On R I establish my working directory by using the navigation tab to the right, once I found the folder where my aligned sequences are I click on the blue icon and choose "Set as working directory", this means that you can call any file inside that working directory just by putting its name.

#Another way to do it is typing this command (on Mac):

> setwd("/Users/admin/Documents/botany563_/Phyloproject/")

#You can always see what files are in your working directory by typing: dir()

#Then used 
dna <- fasta2DNAbin(file="http://adegenet.r-forge.r-project.org/files/usflu.fasta")

#Change "http:......" With the name of your sequence file, in my case is "vanilla_rbcl-aligned.fast"

#So the command looks like:
> dna <- fasta2DNAbin(file="vanilla_rbcl-aligned.fasta")

[March 9]
==============================================
#How to change headers
I used this command to change my sequences headers before the alignment
the idea was to had this headers changed on my allignment, that didnt work

 > new_name<-word(old_name,start=1,end=3,sep=fixed(" "))
 > ref2<-data.frame(old_name,new_name)
> rename.fasta(infile = "vanilla_rbcl.fasta",ref_table = ref2,outfile = "renamed.vanilla_rbcl.fasta")
> clean.fasta.name("renamed.vanilla_rbcl.fasta")
================================================

#msa package
How to run clustalW, ClustalO, muscle 
Follow this link for installation of msa
http://bioconductor.org/packages/release/bioc/vignettes/msa/inst/doc/msa.pdf
#set up your directory first
#load the following libraries

>library(ape)
library(msa)
library(seqinr)

> myseq_rbcl<-readDNAStringSet("vanilla_matk/rbcl.fasta")
> muscle<-msa(myseq_rbcl, "Muscle")
> muscle2<-msaConvert(muscle, type="seqinr::alignment")
> d<-dist.alignment(muscle2)
> muscle_tree<-nj(d)
>plot(muscle_tree)

#I've tried to used msa package with others as ape and phangorn, that didn't work, aparently the conversion feature on msa package does not transformr well DNAbin objects 

######DISTANCE-BASED METHOD
#That last command seems not to work properly so instead I ran this one:

> dna_vanilla_rbcl <-read.dna("vanilla_rbcl-aligned.fasta")
> D <- dist.dna(dna_vanilla_rbcl, model = "TN93")
> tre<- nj(D)
> tre <- ladderize(tre)
> plot(tre, cex=.6)


######PARSIMONY-BASED METHODS

library(ape)
library(adegenet)
library(phangorn)

#load the alligned fasta file
> dna_vanilla_rbcl<-read.dna("vanilla_rbcl-aligned.fasta")
> dna2_vanilla_rbcl<-as.phyDat(dna_vanilla_rbcl)
> tre.ini <- nj(dist.dna(dna_vanilla_rbcl,model="raw"))
> parsimony(tre.ini, dna2_vanilla_rbcl)
> tre.pars <- optim.parsimony(tre.ini, dna2_vanilla_rbcl)
> plot(tre.pars, cex=0.6)

# installing raxml
First download the software from https://github.com/amkozlov/raxml-ng
then place the folder into your "software" folder

in the terminal: 
./raxml-ng -h
git clone https://github.com/amkozlov/ng-tutorial.git

you'll need to be placed inside the raxml folder and then run:
#this is for running the toy data "bad.fa""
./raxml-ng -check -msa ng-tutorial/bad.fa -model GTR+G
./raxml-ng --check --msa ng-tutorial/bad.fa.raxml.reduced.phy --model GTR+G

to my understanding this both commands are checking the MSA, we are using here a toy data "ba.fa" of 6 characters. After the first command it seems it fixed a couple of errors in the alignment, after the second run our alignment is better

for large data we run: 
./raxml-ng --parse --msa ng-tutorial/prim.phy --model GTR+G

now lets infer the ML tree
./raxml-ng --msa ng-tutorial/prim.phy --model GTR+G --prefix T3 --threads 2 --seed 2
./raxml-ng --msa ng-tutorial/prim.phy --model GTR+G --prefix T3-myseed --threads 2 --seed 3162021

im not sure what the second command does, the program generate three new files: best ML tree, all ML trees and optimized model, at least the firs two can be seen in FigTree

UPDATE: according to a google search: different runs with different seeds may or may not result in different trees, this depends on the dataset,
what is guaranteed is that when you specify the same parsimony seed you'll obtain the same starting trees

#running RAXML on my data
im probably doing something wrong but i couldnt run my aligned data outside the raxml folder so i copied my aligned fasta file into the folder and then run:

./raxml-ng -check -msa vanilla_rbcl-aligned.fasta -model GTR+G
./raxml-ng --check --msa vanilla_rbcl-aligned.fasta.raxml.reduced.phy --model GTR+G

./raxml-ng --msa vanilla_rbcl-aligned.fasta.raxml.reduced.phy --model GTR+G --prefix T3 --threads 2 --seed 2

./raxml-ng --msa vanilla_rbcl-aligned.fasta.raxml.reduced.phy --model GTR+G --prefix T3-museed --threads 2 --seed 3162021

#INSTALLING IQ-TREE
  #download the program into your software folder
  http://www.iqtree.org/#download

  then run an example run 
bin/iqtree2 -s example.phy

#running my data using IQ-TREE
bin/iqtree2 -s vanilla_rbcl_aligned.fasta


#MrBAYES

#INSTALLATION (mac)
Download MrBayes (using the terminal)
I followed this tutorial https://github.com/NBISweden/MrBayes/blob/develop/INSTALL
1. Install Homebrew here https://brew.sh/ (it take a while)
2. use the following commands
    brew tap brewsci/bio
    brew install mrbayes
    
#IMPORT DATA
For running MrBayes you need a nexus file, 
use CLUSTALW to export an aligned NEXUS file
run on the terminal: 
clustalw2 -ALIGN -INFILE=vanilla_rbcl.fasta -OUTFILE=VANILLA_RBCL-aligned -OUTPUT=NEXUS

##CREATE A MMBLOCK
Create a txt file and copy the following text, make sure to create the file inside the same folder where your sequences are

begin mrbayes;
 set autoclose=yes;
 prset brlenspr=unconstrained:exp(10.0);
 prset shapepr=exp(1.0);
 prset tratiopr=beta(1.0,1.0);
 prset statefreqpr=dirichlet(1.0,1.0,1.0,1.0);
 lset nst=2 rates=gamma ngammacat=4;
 mcmcp ngen=10000 samplefreq=10 printfreq=100 nruns=1 nchains=3 savebrlens=yes;
 outgroup Anacystis_nidulans;
 mcmc;
 sumt;
end;

now fused both documents: nexus file+mmblock using:
  cat VANILLA_RBCL-aligned mbblock.txt>VANILLA_RBCL-mb.nex

and finally run MrBAYES
  mb VANILLA_RBCL-mb.nex

an error pop up, because i didnt edit the txt file with an OUTGROUP
i decided to include an Epistephium species as an OUTGROUP in my original fasta file and redo all the analysis

Now it worked!
i did the same with my matk sequences #is_alive!






