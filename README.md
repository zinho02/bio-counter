# Paper  
Our paper is the file "Ascending_biological_counters.pdf".  

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
Here we have the two Matlab scripts of a Represssilator implementation. Only the the "repressilatorCounterPaper.m" was used in our paper. The other one could be used to expand this implementation of a counter based on the Repressilator.  
In the end of the file "repressilatorCounterPaper.m" we have 3 plots that were used.  
The first one in line 157: Shows the concentration of each protein in a standard Repressilator implementation.  
The second one in line 160: Shows the States "ladder" as well as the the states 0,1 and 2 described in the paper.  
The third one in line 161: Shows the States "ladder" as well as the the states 3,4 and 5 described in the paper.  
  
## references  
Here we have the main references that we used to have a basis from where to start our paper.  
  
## switch  
Here we have all the switches implementations that we did during the scientific study.  
### initial ideas  
Here we have the Collins toggle switch standard implementation in the file "CollinsToggleSwitch.m".  
The "ourSwitch_1bitCounter.m" was our first try to make a switch.  
### fragile counter  
Here we have better implementations than our initial ideas, but they yet are not stable enough.  
### robust counter  
In the file "Used/ourSwicth_1bitCounter.m" we have our toggle switch with only one input that was used in the paper.  
In the file "Used/3bitCounter.m" we have our implementation of a 3 bit counter used in the paper.    
For all of these models:  
(i)   We have exchanged the name of the repressor ri to pi.  
(ii)  The name of the elements of the second switch have "se" in their beginning.  
(iii) The name of the elements of the third switch have "th" in their beginning.  
