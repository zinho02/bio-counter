# Files

## PulseGenerator
A C++ script to generate a Matlab source code to make a pulse in our switch implementation.
In the file "ent.txt" we have to type the time at each pulse will begin.
The output Matlab code is in the file "saida.txt".
To run this script we just have to run the makefile.

## drawio
Here we have the drawio diagrams that we used and not used in our paper, separeted in the respectives folders.

## images
Here we have the images that we used and not used in our paper, separeted in the respectives folders.

## repressilator
Here we have the two Matlab scripts of a Represssilator implementation. Only the the "repressilatorPlotArtigo.m" was used in our paper. The other one could be used to expand this implementation of a counter based on the Repressilator.
In the end of the file "repressilatorPlotArtigo.m" we have 3 plots that were used.
The first one in line 157: Shows the concentration of each protein in a standard Repressilator implementation.
The second one in line 160: Shows the States "ladder" as well as the the states 0,1 and 2 described in the paper.
The third one in line 161: Shows the States "ladder" as well as the the states 3,4 and 5 described in the paper. 

## sources
Here we have the main sources that we used to have a basis from where to start our paper.

## switch
Here we have all the switches implementations that we did during the scientific study. 
### initial ideas
Here we have the Collins toggle switch standard implementation in the file "switch_pulso_duplo.m".
The "switch_pulso_simples.m" was our first try to make a switch osciilate.
### fragile counter
Here we have better implementations than our initial ideas, but they yet are not stable enough.
### robust counter
