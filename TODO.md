# TODO list

There are more than a few things I should do to make this less stupid.

1.  Currently there are 3 places one might get the 'condition' and 'batch' information
 a. expt$conditions/expt$batches
 b. pData(expt)
 c. expt$design$condition
That is unfortunate and dumb, or more accurately, that is fine, but I need to make clear
which of these is the canonical source when doing analyses and make sure that if I am
going to keep them all, then they need to stay consistent.  Or get rid of the redundant ones.
I am a fan of having a simple slot to see: expt$conditions; but the power of an expressionset
is reduced if one does not use pData(expt), so having it there is good
With that in mind, expt$design is possibly more redunant than worth while.
