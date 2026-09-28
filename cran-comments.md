## Submission

This is a new package.

## R CMD check results

0 errors | 0 warnings | 1 note

* This is a new submission.

## Notes for the reviewers

* The DOIs in the Description field were checked against Crossref.
* 'skedastic' and 'lmtest' are in Suggests: 'skedastic' adds optional
  comparison tests to `run.all.het()`, and 'lmtest' is used only in the tests
  to verify the in-house Goldfeld-Quandt, Breusch-Pagan and White
  implementations.
* The least trimmed squares fits use random subsampling. They are run under a
  fixed seed that is restored afterwards, so results are reproducible and the
  user's random-number stream is not changed.
