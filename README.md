REVISION IN PROGRESS
# 2023_relaxed_selection_project
##Background details
This is the master code for PHENOTYPIC DIVEREGE BETWEEN HATCHERY PINK AND COHO SALMON AND THEIR WILD COUNTERPARTS (in press). It is the culmination of Alex Reich's grad school project work.

## Description of the data and file strucutre

## Sharing/Access information
The data can be accessed via this Dryad link: [ADD LINK HERE] (pending Dryad publication)

# #Code/ Sfotware

_________________________________________
OLD NOTES:
Revisiting code to my masters project. My old project/code are acting up on my new work computer, so creating a new project to clean up and circumnavigate the issue.
I also cleaned up the code quite a bit, but not entirely. It's CLEAN-ER. Also it uses the full coho dataset, not the subsetted coho dataset.
This is a revision of the Relaxed_selection_proj_reproducible_science R project that was my best attempt to make it reproducible when I was in grad school
and drowning in stats exams. This version is cleaner.

This is the readme file.

Masters_thesis_organization.R
- this R file analyzes GSI and egg diameter for the pink 2020, pink 2021, and full coho datasets. It spits out two csv's into the results folder, one for
--GSI and the other for egg. It also spits out the GSI/EGG graph into the Plots folder and results into the results folder.

Bomb_final_simple_dec2023.Rmd
- this rmd analyzes the bomb data. There's some excess exploratory stuff in there, but essentially the file spits out one combined egg energy and lipid 
--csv into the Results folder and a two plots into the Plots folder: a boxplot (for the supplement) and a violin plot (we chose not to use this one, but it exists).

Fig3_lengthcompare_exerything_dec23.R
-prints out the length graph, which, ironically, is now figure 2 in the scientific manuscript. But it was figure 3 in the thesis manuscript.

post_hoc_power_analyses.R
- this was an outline for my post-hoc power analyses. I ended up doing them on Masters-thesis-organization.R for females and Male_morpho_analysis_cleanup_dec_2023.R for males4

Male_morpho_analysis_cleanup_dec_2023.R
- the male linear and geomorph analyses. spits out figures and tables into the Plots and Results folders.
