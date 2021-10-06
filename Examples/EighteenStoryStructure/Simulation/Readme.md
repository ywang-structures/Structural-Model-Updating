"EighteenStoryStructureSim.m" 
* This is the main function to perform FE model updating of the 18-story structure from 100 starting points;

"\Results" 
* This is a folder including the model updating results of the 18-story structure and the code to plot results
* Please run "PostProcessMAC.m","PostProcessEigVec.m", and "PostProcessMDR.m" to plot results
* Followings are the descriptions of each result file:  
     1. EighteenStoryStructureSim_form1_JACoff_Levenberg-Marquardt.mat     : applying L-M algorithm on MAC value formulation using numerical Jacobian 
     2. EighteenStoryStructureSim_form1_JACon_Levenberg-Marquardt.mat      : applying L-M algorithm on MAC value formulation using analytical Jacobian
     3. EighteenStoryStructureSim_form2_JACoff_Trust-region-reflective.mat : applying TRR algorithm on eigenvector difference formulation using numerical Jacobian 
     4. EighteenStoryStructureSim_form2_JACon_Trust-region-reflective.mat  : applying TRR algorithm on eigenvector difference formulation using analytical Jacobian 
     5. EighteenStoryStructureSim_form3_JACon_Levenberg-Marquardt.mat      : applying L-M algorithm on modal dynamic residual formulation using analytical Jacobian 
     6. EighteenStoryStructureSim_form3_JACon_Trust-region-reflective.mat  : applying TRR algorithm on modal dynamic residual formulation using analytical Jacobian
