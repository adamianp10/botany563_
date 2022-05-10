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



====================================================
## MP ANALYSIS

#I decided to used the CIPRES website which have several Phylogenetic programs implemented in their server
https://www.phylo.org/

#First tyou upload the input files as an aligned data matrix, in my case the fasta aligned files of matk, rbcl, nrITS, and then selec the tool, i choose the MPboot XSEDE with the following parameters:

runtime=0.5h
seed_val=12345
sepecyIterations=auto
specifyUboots=1000

#the program does not take so long to procces the imput data and provided several output files including a .mpboot file and a parsimonious tree (.parstree) as well as a concensus tree (.contree). The mpboot file works as a log with have all the running information. For the analysis i decided to use the consensus tree.

## ML ANALYSIS

#Here i decided to go for IQ-tree since is really easy to use and has the plus that incorporate ModelFinder which test among different subtitution models and picks automatically the one with the highest AIC score. In addition has bult-in UFbootstrap. 
#For runnning my data i used the following command (in this case for my rbcl fasta file):

bin/iqtree2 -s vanilla_rbcl_aligned.fasta -bb 1000

##bb refers to the bootstrap iterations which in this case is set up for 1000 repetitions. 

#after running the program, it produced several files including an iqtree which functions a a log with all the parameter of the running including the different models tested as well as a ML and a consensus tree, for the purpose of this project i choose the latter for interpretation. 

## BI ANALYSIS

#Installing Beast/Beauti
Download from here http://www.beast2.org/

#After installing Beast, a pacjage of several tools installed in my computer including BEAST, BEAUTi, LogCombiner, TreeAnnotator and DensiTree. Tracer will need to and can be download from here 
https://github.com/beast-dev/tracer/releases/tag/v1.7.2

#To get a good idea how BEAST worked i followed the beggenigers tutorial here:
https://taming-the-beast.org/tutorials/Introduction-to-BEAST2/#fig:tracer_joint

#I found that BEAST has a bult-in packacgae for model sites called bModelTest that can be installed going to BEAUTi to file/Lauch apps. I decided to used that packacge for my data following this two tutorials:

https://github.com/Taming-the-BEAST/Substitution-model-averaging
https://github.com/BEAST2-Dev/bModelTest/wiki

#In general terms the pipeline goes as this:
- Import your data to BEAUTI (fasta aligned files)
- decided if you want to link your files (this is only needed if you had more than one DNA region, or i you thing certaing regons may have different assumptions (clocl model, site model, etc).
- Choose the site model, the program only give you a couple that according to some FAQ you can modified to the different other models. In my case i installed first the bModelTest which choose the model for me and selected this option in the site Model tab. 
- Then, set priors, which in my case was limited to leave the values as default.
- Set up the MCMC chain lenght, this vary among data, but a good staring point is 10 million. In addition you need to enter the store every, tracelog, screenlog and treelog values which can be 1000 in this example
-Now save this settings as a xml file.
-Then, open BEAST and enter as an output the xml file you generated and run.
- The program will generate several files including a log and a trees files
- Open the log in Tracer and look for the ESS values, anything below 200 will be in yellow or red which means low convergence, which in few words means, that the MCMC needs more time.
- If you are using the bModeltest look for the BMT_ModelIndicator which will give the time the posterior value the chain spend more time, the bar with the highest value correspon to the more support site model. The X axis only gives you numbers, you'll need to refer to the puplication of the package or the tutorial for an interpretation of the model, for example a bar on 17 means K81 model. 
- Then, open TreeAnnotor and load the tree file, the burnin percentage and the PP limit depends on the user, but as a thubm rule you can use a 20% and 0.5 PP, which is the threshold i used on my data. Following the tutorial you'll need to change the node heights to mean hts. 
- Finally this will produce a tree file that you can visualize in FigTree. the posterior values are the ones that the program would ask you to introduce at the begginig and later can be visualize using the branches label tab.

#My first attemps was not succesfull, with ESS values way below 200, which according to the tutorial is not desirable. After doing some research I found that setting up the values of MCMC to 20 and 40 million may results in convergence and higher ESS values, i intended with 20 million which worked for matk adn rbcl but not for nrITS. #After some testing i ended up using 70 million iterations for my nrITS data which gave me high ESS values. This took almost 14 hours to run! 

#TIP LABELS

#I realized that in order to modify the labels i needed to have my  fasta file with the correct labels since the beggining of my pipeline, that means before running the allignment. After several hours looking for an easy way to change my labels i found this small tutorial to do it, Unfortunaly does not pull the posterior or bootstrap values, but worth trying. I endded up changing my tip labels manually.
#source:  https://www.researchgate.net/post/How_to_edit_tip_labels_in_MEGA
#ScripT:
> vanilla_its<-read.tree("final_treev2")
> tree2<-read.nexus("final_ITS_copy.tree")
> tips<-read.csv("tip_labels.csv$species")
> tips<-read.csv("tip_labels.csv")
> newtree<-phylotools::sub.taxa.label(tree2,tips)
> ape::write.tree(newtree, file="Newtree.txt")

#SORRY_FOR_MY_TERRIBLE_GRAMMAR



