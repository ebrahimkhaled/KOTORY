## Submission

This is a new package.

## Test environments

* local: Windows 11, R 4.4.1
* win-builder: R-devel (2026-09-25 r90590 ucrt)
* GitHub Actions: macOS (release), Windows (release), Ubuntu (devel, release, oldrel-1)

## R CMD check results

0 errors | 0 warnings | 1 note

* checking CRAN incoming feasibility ... NOTE
  New submission.
  "Possibly misspelled words in DESCRIPTION": Breusch, Goldfeld, Hartley,
  Hartley's, Quandt and Rousseeuw are the surnames of the authors of the
  cited methods, and are spelled correctly.

Locally there is also "checking for future file timestamps ... NOTE: unable to
verify current time", which comes from the offline check machine.

## Notes for the reviewers

* The DOIs in the Description field were checked against Crossref.
* 'skedastic' and 'lmtest' are in Suggests: 'skedastic' adds optional
  comparison tests to `run.all.het()`, and 'lmtest' is used only in the tests
  to verify the in-house Goldfeld-Quandt, Breusch-Pagan and White
  implementations.
* The least trimmed squares fits use random subsampling. They are run under a
  fixed seed and the previous random-number state is restored afterwards, so
  results are reproducible and the user's random-number stream is not changed.
  Some functions of 'skedastic' call set.seed() internally; `run.all.het()`
  restores the user's state after calling them as well.
