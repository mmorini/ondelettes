# ondelettes
Ondelettes project python scripts

1. PREPROCESSING

   1a. Parse WoS "savedrecs" data to obtain ".dat" files
   1b. Handclean files (Astrophysics Ade PAR articles*)

   1c. Make temporal slices

2. COMMUNITY DETECTION

   2. Run Louvain algorithm (atoms version) on each slice

3. OBTAIN IdCards (optional)

   3. Make pdf files out of tex files

----------------------------------------------------------
1.-3.: bash waves.sh [check directories]
----------------------------------------------------------


4. BUILD HISTORICAL NETWORK

   4. Link temporal slices (jaccard distances, colors, quality metric)
      we see more links than usual, workingon ATOMS, not usual comms

----------------------------------------------------------
4.: python jaccard_history.py [check directories]
----------------------------------------------------------

5. [does not exist anymore]

6. add_attributes_colors_positions, final cleaning 

----------------------------------------------------------
6. finalize_graph.py [check directories]
----------------------------------------------------------
