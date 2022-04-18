# OpenFASTWrappers

A repository containing julia wrappers for standalone openFAST modules.

## Quick Start for Developing Julia Packages

- include("./path/modulename.jl/src/modulename.jl") is like a copy and paste and will reload the module every time the script is run

- installing a module by navigating to the module's working directory (and in an interactive julia session), running "] dev ." will install it in dev mode and will pick up changes made in the local package location every time julia restarts.

- installing normally ("] add url/or/path/2/package") installs the package in the "~/.julia" folder and requires the url/path original location to be updated, and a "] update mypackage" to get it to pick up changes.

- Each package has its own "environment" like conda environment.  To activate that environment and add dependencies, navigate to the working directory and run "] activate ." Then, when you "] add dependencyname" it will add it to the Project.toml and Manifest.toml automatically.  To switch back to the general environment, run "] activate"

- Add code in the src/ folder with an include("./yournewfile.jl") in the OpenFASTWrapers.jl file

- Add tests and test data in the test/ folder, add include statements for your tests in the runtests.jl file.
