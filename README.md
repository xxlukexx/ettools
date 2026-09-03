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

## Normalized gaze coordinates

`etGazeData` treats a sample from an eye as missing when either normalized
coordinate is non-finite or outside the closed interval `[0, 1]`. By default,
both coordinates for that eye are stored as `NaN`, and binocular averages use
only eyes that remain valid. Explicit missing-data masks are also honoured.

For forensic access to source coordinates, construct an object with
`'PreserveInvalidCoordinates', true`. This retains the original per-eye values
while still marking them as missing and excluding them from binocular averages.
