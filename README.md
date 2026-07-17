# ettools

Shared MATLAB eye-tracking utilities used by Task Engine analyses.

`etPreprocess.m` is the established ECK/Task Engine main-buffer preprocessor.
It was recovered into this canonical repository from
`gen2020/matlab/ECKAnalyse` on 17 July 2026 after the PIP conformance run
demonstrated that `ms` otherwise depended on Gen2020 being manually added to
the MATLAB path. The associated time-buffer helpers used by `etResample2`
were recovered at the same time. Study runners must use the canonical
workspace path setup; they must not add a study checkout to satisfy this
dependency.
