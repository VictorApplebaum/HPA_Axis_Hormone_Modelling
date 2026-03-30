# HPA_Axis_Hormone_Modelling
Contains code for the currently unpublished statistical modelling framework developed for my PhD. We examine the relationship between 7 HPA axis hormones using a Bayesian hierarchical multiplicative vector autoregressive (mVAR) model.

Instructions:
1. Run DataPreprocessing.R, which does all the processing, saving df_full.RData which will be used in future. As written, it currently is set up to run 20 folds, however we use 10 for the results in the thesis/paper
2. Run MainRun.R to get the main results for the project, including the heatplot, using the model described in Model.stan
3. Run FoldTraining.R to run the folds, numbered using FoldNo at the top of the document. I did these in parallel with the following bash code:
for i in {1..20}
do
  Rscript FoldTraining.R --FoldNo=$i &
done
wait
4. Run FoldTesting.R similarly to make test predictions
5. For evaluation of results from the FoldTesting.R runs, use FoldEvaluating.R
6. To get plots of predictions made for individuals by specific folds on the tasks, use PredictionPlotter.R
7. To get results for reduced priors (sensitivity analysis in the supplementary) use MainRunReducedLambdaPriors.R.
