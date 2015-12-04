Dear Robert,

I regret that we were not able to address your concerns directly.

However, I would like to point out the following:
  - It is not necessary to propagate down the GO terms to leaf terms. (you just need to change gene*GO matrix to gene*GO-of-your-choice)
- The purpose of using leaf terms is to try to make the terms used be as independent as possible, which is a better practice in machine learning aspects.
- We had tested on original GO terms  (not propagated down) and found that using leaf nodes performed better.
+ We have not tested on GOSlim. Thank you for the suggestion.
- We are not saying that GO term propagation (up or down) is generally correct way to use GO term. However, assumption is that the significant terms will be enriched (however, results may have false positives, which is unavoidable in any methods)
+ This cannot be tested based on one gene (we had removed samples with small number of mutated genes) :
  + I realize we should have stated this more clearly in the paper.
