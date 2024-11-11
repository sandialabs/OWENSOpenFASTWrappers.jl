# OWENSOpenFASTWrappers

[![](https://img.shields.io/badge/docs-stable-blue.svg)](https://sandialabs.github.io/OWENSOpenFASTWrappers.jl)
[![](https://img.shields.io/badge/docs-dev-blue.svg)](https://sandialabs.github.io/OWENSOpenFASTWrappers.jl/dev)
![](https://github.com/sandialabs/OWENSOpenFASTWrappers.jl/workflows/CI/badge.svg)
[![codecov](https://codecov.io/gh/sandialabs/OWENSOpenFASTWrappers.jl/graph/badge.svg?token=1egaz9f5cT)](https://codecov.io/gh/sandialabs/OWENSOpenFASTWrappers.jl)

* *Documentation:* [sandialabs.github.io/OWENSOpenFASTWrappers.jl/](https://sandialabs.github.io/OWENSOpenFASTWrappers.jl/dev)
* *Code:* [github.com/sandialabs/OWENSOpenFASTWrappers.jl](https://github.com/sandialabs/OWENSOpenFASTWrappers.jl)

This is part of the Sandia National Labs OWENS toolkit (Onshore/Offshore Wind/Water Energy Simulator).  This package provides a frontend julia wrapper to many of the popular OpenFAST libraries, including AeroDyn, HydroDyn, MoorDyn, InflowWind, and Turbsim.  The package can be used standalone, as shown in the test cases, or in conjunction with OWENS.jl.  The OWENSOpenFAST_jll.jl package was also created, which this is linked to, which provides cross compiled binaries for the major operating systems and significantly improves the installation experience.  Please also note the API reference in these docs, which gives a searchable index of all of the functions with their inputs and outputs. Please make any pull requests against the dev branch.

## Installation
]add git@github.com:sandialabs/OWENSOpenFAST_jll.jl.git
]add git@github.com:sandialabs/OWENSOpenFASTWrappers.jl.git