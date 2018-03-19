This folder contains R scripts to run analysis for PHENOTYPES and COMMUNITY & ECOSYSTEMS. Results have to be saved in a *results* folder. Then, we used R packages *rstanarm* and *bayesplot* (see R script *PosteriorCheck*) to perform analysis, predictions and diagnostics.

Models comparison were performed using the R package *loo*.

Note that the R script "ANALYSIS_COMMUNITY_ECOSYSTEM.R" do not replicate models for all datasets (i.e, With fish, without fish, Surber,...). Instead, one have to select the dataset to analyse (see lines 35-38).

## PHENOTYPES

- Exp_1_Mass_Growth_Hatchery_Conditions => Fig 2a and 2b + Extended data Table 1
Response variables: Body_Mass and Growth_Rate_SGRW

- Exp_1_Mass_Growth_Experimental_Streams => Fig 2a and 2b + Extended data Table 1
Response variables: Body_Mass, Growth_Rate_SGRW_Before_Release and Growth_Rate_SGRW_After_Release

- Exp_1_Morphology_Experimental_Streams => Fig 2c + Extended data Table 1 + Extended data Figure 2
Response variables: WARP1 and WARP2

- Exp_2_Activity_Open_Field_Tests => Fig 2d + Extended data Table 1 
Response variables: Activity

- Exp_2_Movement_Habitat_Use_Stream_Mesocosms => Fig 2e + Extended data Table 1
Response variables: Movement and Habitat_Use

- Exp_3_Excretion_Hatchery_Conditions => Extended data Table 1
Response variables: N_Excretion and P_Excretion 

- Exp_3_Excretion_Stream_Mesocosms => Fig 2f + Extended data Table 1
Response variables: N_Excretion and P_Excretion 


## COMMUNITY AND ECOSYSTEMS

- Exp_1_Community_Ecosystem_With_Fish_Experimental_Streams => Fig 3, Extended data Figure 4 + Extended data Table 2
Response variables: Rhyacophilidae, Polycentropodidae, Total_Prim_Cons, Chironomidae, Hydropsychidae, Baetidae, Simuliidae and Planorbidae, Primary_production, Total_Decomposition and Microbial_Decomposition 

- Exp_1_Surber_With_Fish_Experimental_Streams => Extended data Figure 3+ Extended data Table 3
Response variables:Rhyacophilidae, Polycentropodidae, Total_Prim_Cons, Chironomidae, Hydropsychidae, Baetidae, Simuliidae and Planorbidae 


- Exp_1_Community_Ecosystem_Without_Fish_Experimental_Streams => Extended data Figure 5 + Extended data Table 2
Response variables: Rhyacophilidae, Polycentropodidae, Total_Prim_Cons, Chironomidae, Hydropsychidae, Baetidae, Simuliidae and Planorbidae, Primary_production, Total_Decomposition and Microbial_Decomposition 

- Exp_1_Surber_Without_Fish_Experimental_Streams => Extended data Table 3
Response variables: Rhyacophilidae, Polycentropodidae, Total_Prim_Cons, Chironomidae, Hydropsychidae, Baetidae, Simuliidae and Planorbidae 

- Exp_3_Community_Ecosystem_Stream_Mesocosms => Table 1, Extended data Figure 6 and Extended data Table 4
Response variables: Chironomidae, Primary_Production, Total_Decomposition and Microbial_Decomposition 

- Exp_3_Excretion_Population_Stream_Mesocosms => Table 1 and Extended data Table 4
Response variables: N_Excretion and P_Excretion

- Exp_3_Filamentous_Algae_Stream_Mesocosms => Table 1 and Extended data Table 4
Response variables: Filamentous_Algae

- Exp_3_Nutrient_Fluxes_Stream_Mesocosms=> Table 1 and Extended data Table 4
Response variables: NH4_concentration and PO4_concentration

